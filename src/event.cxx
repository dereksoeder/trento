// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <algorithm>
#include <cmath>
#include <complex>

#include <boost/program_options/variables_map.hpp>

#include "nucleus.h"

namespace trento {

namespace {

typedef std::complex<double> complex_t;

constexpr double TINY = 1e-12;

// Generalized mean for p > 0.
// M_p(a, b) = (1/2*(a^p + b^p))^(1/p)
inline double positive_pmean(double p, double a, double b) {
  return std::pow(.5*(std::pow(a, p) + std::pow(b, p)), 1./p);
}

// Generalized mean for p < 0.
// Same as the positive version, except prevents division by zero.
inline double negative_pmean(double p, double a, double b) {
  if (a < TINY || b < TINY)
    return 0.;
  return positive_pmean(p, a, b);
}

// Generalized mean for p == 0.
inline double geometric_mean(double a, double b) {
  return std::sqrt(a*b);
}

inline int get_max_eccentricity_m(const std::vector<EventQuantity>& quantities) {
  int max_m = -1;
  bool foundmn = false;

  for (auto qty : quantities) {
    if (EventQuantity_GetClass(qty) == EventEpsilon_mn) {
      int m = EventQuantity_GetSubscript1(qty);
      int n = EventQuantity_GetSubscript2(qty);

      max_m = std::max(max_m, m);
      foundmn = (foundmn || (m != n));
    }
  }

  return (foundmn ? max_m : -1);
}

inline int get_max_eccentricity_n(const std::vector<EventQuantity>& quantities) {
  int max_n = -1;

  for (auto qty : quantities) {
    if (EventQuantity_GetClass(qty) == EventEpsilon_mn) {
      int n = EventQuantity_GetSubscript2(qty);
      max_n = std::max(max_n, n);
    }
  }

  return max_n;
}

}  // unnamed namespace

Event::Event(const VarMap& var_map)                                    // this constructor assumes that only quantities being displayed need to
  : Event(var_map, var_map["columns"].as<EventQuantityList>().values)  //   be computed; retained for code compatibility but should not be used
{ }

// Determine the grid parameters like so:
//   1. Read and set step size from the configuration.
//   2. Read grid max from the config, then set the number of steps as
//      nsteps = ceil(2*max/step).
//   3. Set the actual grid max as max = nsteps*step/2.  Hence if the step size
//      does not evenly divide the config max, the actual max will be marginally
//      larger (by at most one step size).
Event::Event(
  const VarMap& var_map,
  const std::vector<EventQuantity>& required_quantities)
    : norm_(var_map["normalization"].as<double>()),
      dxy_(var_map["grid-step"].as<double>()),
      nsteps_(std::ceil(2.*var_map["grid-max"].as<double>()/dxy_)),
      xymax_(.5*nsteps_*dxy_),
      TA_(boost::extents[nsteps_][nsteps_]),
      TB_(boost::extents[nsteps_][nsteps_]),
      TR_(boost::extents[nsteps_][nsteps_]),
      max_eccentricity_m_(get_max_eccentricity_m(required_quantities)),
      max_eccentricity_n_(get_max_eccentricity_n(required_quantities)) {
  // Choose which version of the generalized mean to use based on the
  // configuration.  The possibilities are defined above.  See the header for
  // more information.
  auto p = var_map["reduced-thickness"].as<double>();

  if (std::fabs(p) < TINY) {
    compute_reduced_thickness_ = [this]() {
      compute_reduced_thickness(geometric_mean);
    };
  } else if (p > 0.) {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return positive_pmean(p, a, b); });
    };
  } else {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return negative_pmean(p, a, b); });
    };
  }
}

void Event::compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
                    const NucleonCommon& nucleon_common) {
  // Reset npart; compute_nuclear_thickness() increments it.
  npart_ = 0;
  compute_nuclear_thickness(nucleusA, nucleon_common, TA_);
  compute_nuclear_thickness(nucleusB, nucleon_common, TB_);
  compute_reduced_thickness_();
  compute_observables();

  failuresA_ = nucleusA.failures();
  failuresB_ = nucleusB.failures();
}

namespace {

// Limit a value to a range.
// Used below to constrain grid indices.
template <typename T>
inline const T& clip(const T& value, const T& min, const T& max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

}  // unnamed namespace

void Event::compute_nuclear_thickness(
    const Nucleus& nucleus, const NucleonCommon& nucleon_common, Grid& TX) {
  // Construct the thickness grid by looping over participants and adding each
  // to a small subgrid within its radius.  Compared to the other possibility
  // (grid cells as the outer loop and participants as the inner loop), this
  // reduces the number of required distance-squared calculations by a factor of
  // ~20 (depending on the nucleon size).  The Event unit test verifies that the
  // two methods agree.

  // Wipe grid with zeros.
  std::fill(TX.origin(), TX.origin() + TX.num_elements(), 0.);

  // Deposit each participant onto the grid.
  for (const auto& nucleon : nucleus) {
    if (!nucleon.is_participant())
      continue;

    ++npart_;

    // Get nucleon subgrid boundary {xmin, xmax, ymin, ymax}.
    const auto boundary = nucleon_common.boundary(nucleon);

    // Determine min & max indices of nucleon subgrid.
    int ixmin = clip(static_cast<int>((boundary[0]+xymax_)/dxy_), 0, nsteps_-1);
    int ixmax = clip(static_cast<int>((boundary[1]+xymax_)/dxy_), 0, nsteps_-1);
    int iymin = clip(static_cast<int>((boundary[2]+xymax_)/dxy_), 0, nsteps_-1);
    int iymax = clip(static_cast<int>((boundary[3]+xymax_)/dxy_), 0, nsteps_-1);

    // Add profile to grid.
    for (auto iy = iymin; iy <= iymax; ++iy) {
      for (auto ix = ixmin; ix <= ixmax; ++ix) {
        TX[iy][ix] += nucleon_common.thickness(
          nucleon, (ix+.5)*dxy_ - xymax_, (iy+.5)*dxy_ - xymax_
        );
      }
    }
  }
}

template <typename GenMean>
void Event::compute_reduced_thickness(GenMean gen_mean) {
  double sum = 0.;
  double ixcm = 0.;
  double iycm = 0.;

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      auto t = norm_ * gen_mean(TA_[iy][ix], TB_[iy][ix]);
      TR_[iy][ix] = t;
      sum += t;
      // Center of mass grid indices.
      // No need to multiply by dxy since it would be canceled later.
      ixcm += t * static_cast<double>(ix);
      iycm += t * static_cast<double>(iy);
    }
  }

  multiplicity_ = dxy_ * dxy_ * sum;
  ixcm_ = ixcm / sum;
  iycm_ = iycm / sum;
}

void Event::compute_en(int max_n) {
  // Compute eccentricity.
  if (max_n < MinEccentricityHarmonic || max_n > MaxEccentricityHarmonic)
    throw std::invalid_argument{"max_n"};

  // Simple helper class for use in the following loop.
  struct EccentricityAccumulator {
    complex_t z = 0.;
    double wt = 0.;  // weight
    double finish() const  // compute final eccentricity
    { return std::abs(z) / std::fmax(wt, TINY); }
  } en[NumEccentricityHarmonics];

  for (int iy = 0; iy < nsteps_; ++iy) {
    auto y = static_cast<double>(iy) - iycm_;
    auto y2 = y*y;

    for (int ix = 0; ix < nsteps_; ++ix) {
      const auto& t = TR_[iy][ix];
      if (t < TINY)
        continue;

      // Compute `r` relative to the CM.
      auto x = static_cast<double>(ix) - ixcm_;
      auto r = std::sqrt(x*x + y2);

      // The eccentricity harmonics are weighted averages of r^n*exp(i*n*phi)
      // over the entropy profile (reduced thickness).  Note that:
      //
      //   r^n * exp(i*n*phi)  =  r^n * exp(i*phi)^n
      //                       =  (r*exp(i*phi))^n
      //                       =  (x + i*y)^n
      //
      // The Event unit test verifies that this method agrees with the
      // trigonometric function method.
      auto z    = complex_t{x, y};
      auto enz  = t * z;
      auto enwt = t * r;

      for (int i = MinEccentricityHarmonic; i <= max_n; ++i) {
        en[i - MinEccentricityHarmonic].z  += (enz  *= z);
        en[i - MinEccentricityHarmonic].wt += (enwt *= r);
      }
    }
  }

  for (int n = MinEccentricityHarmonic; n <= max_n; ++n)
    eccentricity_[n] = eccentricity_mn_[std::make_pair(n,n)] = en[n - MinEccentricityHarmonic].finish();
}

void Event::compute_emn(int max_m, int max_n) {
  // Compute eccentricity.
  if (max_m < MinEccentricityHarmonic || max_m > MaxEccentricityHarmonic)
    throw std::invalid_argument{"max_m"};

  if (max_n < MinEccentricityHarmonic || max_n > MaxEccentricityHarmonic)
    throw std::invalid_argument{"max_n"};

  // Simple helper class for use in the following loop.
  complex_t emn[NumEccentricityHarmonics][NumEccentricityHarmonics] = { };  // important: initializes to zeroes
  double weights[NumEccentricityHarmonics] = { };

  for (int iy = 0; iy < nsteps_; ++iy) {
    auto y = static_cast<double>(iy) - iycm_;
    auto y2 = y*y;

    for (int ix = 0; ix < nsteps_; ++ix) {
      const auto& t = TR_[iy][ix];
      if (t < TINY)
        continue;

      // Compute `r` relative to the CM.
      auto x = static_cast<double>(ix) - ixcm_;
      auto r = std::sqrt(x*x + y2);

      auto wt = t * r;  // initialize to `(t * r)` and multiply `r` again at top of `m` loop, since `m` starts at 2
      for (int m = MinEccentricityHarmonic; m <= max_m; ++m)       // (note: `MinEccentricityHarmonic = 2` assumed)
        weights[m - MinEccentricityHarmonic] += (wt *= r);

      // The eccentricity harmonics are weighted averages of `r^m*exp(i*n*phi)`
      // over the entropy profile (reduced thickness).  Note that:
      //
      //   r^m * exp(i*n*phi)  =  r^m * exp(i*phi)^n
      //                       =  r^(m-n) * (r*exp(i*phi))^n
      //                       =  r^(m-n) * (x + i*y)^n
      //
      auto z = complex_t{x, y};
      auto zfactor = t * z;  // initialize to `(t * z)` and multiply `z` again at top of `n` loop, since `n` starts at 2

      for (int n = MinEccentricityHarmonic; n <= MaxEccentricityHarmonic; n++) {  // (note: `MinEccentricityHarmonic = 2` assumed)
        zfactor *= z;

        double rfactor;

        if (r > TINY) {  // ignore this pixel if `r ~ 0` here; it would be nugatory in an integral but causes problems for our sum
          rfactor = 1.;  // initialize to 1 and divide `r` at top of `m` loop, since `(m - n)` starts at -1 and counts down

          for (int m = n - 1; m >= MinEccentricityHarmonic; --m) {
            rfactor /= r;
            emn[m - MinEccentricityHarmonic][n - MinEccentricityHarmonic] += rfactor * zfactor;
          }

          rfactor = 1.;  // now initialize to 1 and multiply `r` at bottom of `m` loop, since `(m - n)` now starts at 0 and counts up

          for (int m = n; m <= MaxEccentricityHarmonic; ++m) {
            emn[m - MinEccentricityHarmonic][n - MinEccentricityHarmonic] += rfactor * zfactor;
            rfactor *= r;
          }
        }
      }
    }
  }

  for (int m = MinEccentricityHarmonic; m <= max_m; ++m) {
    double wt = std::fmax(weights[m - MinEccentricityHarmonic], TINY);

    for (int n = MinEccentricityHarmonic; n <= max_n; ++n)
      eccentricity_mn_[std::make_pair(m,n)] = std::abs(emn[m - MinEccentricityHarmonic][n - MinEccentricityHarmonic]) / wt;

    if (m <= max_n)  // also update `eccentricity_` for `m == n` if applicable
      eccentricity_[m] = std::abs(emn[m - MinEccentricityHarmonic][m - MinEccentricityHarmonic]) / wt;
  }
}

void Event::compute_observables() {
  if (max_eccentricity_n_ >= 0) {
    if (max_eccentricity_m_ >= 0)
      compute_emn(max_eccentricity_m_, max_eccentricity_n_);
    else compute_en(max_eccentricity_n_);
  }
}

}  // namespace trento
