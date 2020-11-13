// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015-2020 Jonah E. Bernhard, J. Scott Moreland, Derek Soeder
// MIT License

#ifndef EVENTQTY_H
#define EVENTQTY_H

#include <map>
#include <string>
#include <vector>
#include <boost/any.hpp>

#include "fwd_decl.h"

namespace trento {

constexpr int MinEccentricityHarmonic = 2;  // DO NOT CHANGE: implementations in event.cxx assume this is 2
constexpr int MaxEccentricityHarmonic = 8;
constexpr int NumEccentricityHarmonics = (MaxEccentricityHarmonic - MinEccentricityHarmonic + 1);

/// Enumeration of per-event quantities available for output.
typedef enum {
  /// Default value; does not correspond to any quantity.
  EventQuantity_None = 0,

  /// Which event this is in the sequence (0 = first, 1 = second, ...).
  EventNumber,

  /// Impact parameter 'b' of the collision this event represents.
  EventImpactParameter,

  /// Number of nucleon participants, `npart`.
  EventNumParticipants,

  /// Binary collision number, `ncoll`.
  EventNumCollisions,

  /// Integrated reduced thickness, `mult`.
  EventMultiplicity,

  /// Number of nucleon placement failures while constructing the first nucleus.
  EventFailuresA,

  /// Number of nucleon placement failures while constructing the second nucleus.
  EventFailuresB,

  /// x coordinate of center of mass of reduced-thickness grid `TR_`, in units of grid spacing.
  EventGridCMx,

  /// y coordinate of center of mass of reduced-thickness grid `TR_`, in units of grid spacing.
  EventGridCMy,

  /// Base value for all generalized eccentricity harmonics, $\epsilon_{m,n}$ (arXiv:1111.6538).
  EventEpsilon_mn,  // for a specific `m` and `n`, compute `EventEpsilon_mn + (m << 8) + (n << 16)`, or use `EventQuantity_Epsilon_mn`

  /// Base value for arguments $\Phi_{m,n}$ associated with generalized eccentricity harmonics $\epsilon_{m,n}$.
  EventEpsilonArg_mn,  // for a specific `m` and `n`, compute `EventEpsilonArg_mn + (m << 8) + (n << 16)`, or use `EventQuantity_EpsilonArg_mn`

  //
  //---- insert new quantities above this comment, before `EventQuantity_Mask` ----
  //
  EventQuantity_Mask = 0xFF  // use `(qty & EventQuantity_Mask)`, or `EventQuantity_GetClass`, to get the general quantity; specifics (e.g., epsilon subscripts) are in higher bits

} EventQuantity;

static constexpr EventQuantity EventQuantity_Epsilon_mn(int m, int n)    { return static_cast<EventQuantity>(EventEpsilon_mn + (m << 8) + (n << 16)); }
static constexpr EventQuantity EventQuantity_Epsilon_n(int n)            { return EventQuantity_Epsilon_mn(n, n); }

static constexpr EventQuantity EventQuantity_EpsilonArg_mn(int m, int n) { return static_cast<EventQuantity>(EventEpsilonArg_mn + (m << 8) + (n << 16)); }
static constexpr EventQuantity EventQuantity_EpsilonArg_n(int n)         { return EventQuantity_EpsilonArg_mn(n, n); }

static constexpr EventQuantity EventQuantity_GetClass(EventQuantity id) { return static_cast<EventQuantity>(id & EventQuantity_Mask); }
static constexpr int EventQuantity_GetSubscript1(EventQuantity id)      { return (id >>  8) & 0xFF; }
static constexpr int EventQuantity_GetSubscript2(EventQuantity id)      { return (id >> 16) & 0xFF; }

/// Maps numerical identifiers to textual names for event quantities.
extern std::map<EventQuantity, std::vector<std::string>> EventQuantityLabels;

struct EventQuantityList {
  std::vector<EventQuantity> values;

  EventQuantityList() = default;

  EventQuantityList(std::initializer_list<EventQuantity> ids) :
    values(ids)
  { }

  EventQuantityList(const std::vector<EventQuantity>& ids) :
    values(ids)
  { }

  void append(const std::string& token);

  EventQuantityList operator+(EventQuantity id) const {
    EventQuantityList ret(values);
    ret.values.push_back(id);
    return ret;
  }

  inline static bool is_integer(EventQuantity id) {
    // TODO: improve this
    return (id == EventNumber || id == EventNumParticipants || id == EventNumCollisions || id == EventFailuresA || id == EventFailuresB);
  }

  static EventQuantity lookup(const std::string& label);

  static const std::string& lookup(EventQuantity id);

  static EventQuantityList create(const std::string& token) {
    EventQuantityList ret;
    ret.append(token);
    return ret;
  }
} ;

/// Default quantities displayed in Trento output.
extern const EventQuantityList DefaultEventQuantityList;

/// Custom `boost::program_options` validator to pre-parse string argments into `EventQuantityList`s.
void validate(boost::any& v, const std::vector<std::string>& tokens, EventQuantityList *, int);

}  // namespace trento

#endif  // EVENTQTY_H
