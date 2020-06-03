// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015-2020 Jonah E. Bernhard, J. Scott Moreland, Derek Soeder
// MIT License

#include <string>
#include <stdexcept>
#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>

#include "eventbin.h"
#include "eventqty.h"

namespace trento {

namespace {

static void validate_edges(const std::vector<double>& edges) {
  if (edges.size() < 2)
    throw std::invalid_argument{"cannot bin with fewer than two edges"};

  auto it = std::begin(edges);
  double prev = *it;
  double edge = *(++it);
  if (prev == edge)
    throw std::invalid_argument{"bin edges cannot be equal"};
  bool ascending = (prev < edge);

  while (++it != std::end(edges)) {
    prev = edge;
    edge = *it;
    if (prev == edge)
      throw std::invalid_argument{"bin edges cannot be equal"};
    else if ((prev < edge) != ascending)
      throw std::invalid_argument{"bin edges must be monotonic"};
  }
}

static size_t binary_search_ascending_for_discard(const std::vector<size_t>& permwithflag, const std::vector<double>& values, double target) {
  size_t i = 0;
  for (size_t j = permwithflag.size(); i < j; ) {  // keep `j` past the end of the search range to avoid computing `(j - 1)`, which could underflow
    auto k = (i + j - 1) / 2;
    auto value = values[permwithflag[k] >> 1];  // LSB is deletion flag
    if (target > value)
      i = k + 1;
    else if (target < value)
      j = k;
    else return k;
  }

  return i;
}

static void remove_for_discard(const std::vector<size_t>& permwithflag, std::vector<double>& values) {
  size_t n = permwithflag.size();

  size_t i;
  for (i = 0; i < n; ++i)
    if (permwithflag[i] & 1) break;  // LSB is deletion flag

  size_t j = i;
  for ( ; i < n; ++i)
    if (!(permwithflag[i] & 1))
      values[j++] = values[i];

  values.resize(j);
}

static int binary_search_bins_ascending(const std::vector<double>& edges, double target) {
  size_t i = 0;
  for (size_t j = edges.size() - 1; (i + 1) < j; ) {  // ensures `i <= edges.size() - 2`
    auto k = (i + j) / 2;
    auto edge = edges[k];
    if (target < edge)  // smaller (left) edge is inclusive, larger (right) edge is exclusive
      j = k;
    else i = k;
  }

  return ((target >= edges[i] && target < edges[i + 1]) ? static_cast<int>(i) : -1);
}

static int binary_search_bins_descending(const std::vector<double>& edges, double target) {
  size_t i = 0;
  for (size_t j = edges.size() - 1; (i + 1) < j; ) {  // ensures `i <= edges.size() - 2`
    auto k = (i + j) / 2;
    auto edge = edges[k];
    if (target >= edge)  // smaller (right) edge is inclusive, larger (left) edge is exclusive
      j = k;
    else i = k;
  }

  return ((target < edges[i] && target >= edges[i + 1]) ? static_cast<int>(i) : -1);
}

}  // unnamed namespace

EventBinner::EventBinner(const std::vector<EventQuantity>& columns)
  : columns_(columns)
{
  for (auto qty : columns_)
    data_[qty];  // creates entry if it does not exist
}

void EventBinner::set_subsample(size_t count) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  subsample_ = count;
}

void EventBinner::add_axis(EventQuantity qty, std::vector<double>&& edges) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  validate_edges(edges);

  data_[qty];  // creates entry if it does not exist

  // BinAxis::{quantity, edges, ascending, is_percentile} (via constructor)
  axes_.emplace_back( qty, std::forward<std::vector<double>>(edges), (edges.front() < edges.back()), false );
}

void EventBinner::add_axis_percentile(EventQuantity qty, std::vector<double>&& edges) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  validate_edges(edges);

  data_[qty];  // creates entry if it does not exist

  // BinAxis::{quantity, edges, ascending, is_percentile} (via constructor)
  axes_.emplace_back( qty, std::forward<std::vector<double>>(edges), (edges.front() < edges.back()), true );
  delay_binning_ = true;
}

void EventBinner::add_axis_percentile(EventQuantity qty, std::vector<double>&& edges, double include_value, bool in_last_bin) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  validate_edges(edges);

  data_[qty];  // creates entry if it does not exist

  // BinAxis::{quantity, edges, ascending, first_inclusion, last_inclusion} (via constructor)
  axes_.emplace_back( qty, std::forward<std::vector<double>>(edges), (edges.front() < edges.back()), include_value, in_last_bin );
  delay_binning_ = true;
}

void EventBinner::add_axis_percentile(EventQuantity qty, std::vector<double>&& edges, double include_first, double include_last) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  validate_edges(edges);

  data_[qty];  // creates entry if it does not exist

  // BinAxis::{quantity, edges, ascending, first_inclusion, last_inclusion} (via constructor)
  axes_.emplace_back( qty, std::move(edges), (edges.front() < edges.back()), include_first, include_last );
  delay_binning_ = true;
}

void EventBinner::add_discard(EventQuantity qty, bool greater_than, double value) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  data_[qty];  // creates entry if it does not exist

  // DiscardCriterion::{quantity, greater_than, value, is_percentile} (via constructor)
  discard_criteria_.emplace_back( qty, greater_than, value, false );
}

void EventBinner::add_discard_percentile(EventQuantity qty, bool greater_than, double percentile) {
  if (numevents_ > 0)
    throw std::logic_error{"cannot reconfigure binner once processing has started"};

  data_[qty];  // creates entry if it does not exist

  // DiscardCriterion::{quantity, greater_than, value, is_percentile} (via constructor)
  discard_criteria_.emplace_back( qty, greater_than, percentile, true );
  delay_binning_ = true;
}

void EventBinner::start() {
  if (axes_.size() == 0)
    throw std::invalid_argument{"no bin edges were supplied"};

  size_t limit = ((std::numeric_limits<size_t>::max() / 2) / sizeof(decltype(bin_sums_)::value_type)) / std::max(columns_.size(), static_cast<size_t>(1));
  size_t numbins = 1;

  for (const auto& axis : axes_) {
    auto axisbins = (axis.edges.size() - 1);

    if (axisbins > limit)
      throw std::invalid_argument{"too many bins"};

    numbins *= axisbins;
    limit   /= axisbins;
  }

  bin_counts_.resize(numbins);
  bin_sums_.resize(numbins * columns_.size());
}

void EventBinner::process(const Event& event) {
  ++numevents_;

  if (delay_binning_) {
    for (auto& field : data_)
      field.second.push_back(event.get<double>(field.first));

    if (++numevents_pending_ == subsample_ && subsample_ != 0)
      resolve_percentiles();  // clears `delay_binning_`
  }
  else if (!should_discard(event)) {
    std::function<double(EventQuantity)> fn = std::bind(&Event::get<double>, &event, std::placeholders::_1);
    bin_event(fn);
  }
}

void EventBinner::finish() {
  if (delay_binning_) {
    subsample_ = 0;
    resolve_percentiles();
  }
}

bool EventBinner::should_discard(const Event& event) const {
  // all percentiles should have been resolved by this point, so 'cut.is_percentile' should always be false
  for (const auto& cut : discard_criteria_) {
    double value = event.get<double>(cut.quantity);
    if (cut.greater_than ? (value >= cut.value) : (value <= cut.value))  // always discard if value is exactly at cutoff
      return true;
  }

  return false;
}

void EventBinner::bin_event(std::function<double(EventQuantity)> getter) {
  size_t index = 0;
  for (const auto& axis : axes_) {
    index *= (axis.edges.size() - 1);

    auto value = getter(axis.quantity);
    int b = (axis.ascending ? binary_search_bins_ascending : binary_search_bins_descending)(axis.edges, value);
    if (b < 0) return;

    index += static_cast<size_t>(b);
  }

  bin_counts_[index]++;

  index *= columns_.size();
  for (auto qty : columns_) {
    auto value = getter(qty);
    std::pair<double, double>& sums = bin_sums_[index++];
    sums.first  += value;
    sums.second += (value * value);
  }
}

void EventBinner::resolve_percentiles() {
  if (!delay_binning_)
    throw std::logic_error{"no unresolved percentiles"};

  std::vector<size_t> perm;
  size_t nevents = numevents_pending_;

  // apply discard cuts

  if (!discard_criteria_.empty()) {
    std::sort(
      std::begin(discard_criteria_), std::end(discard_criteria_),
      [](const DiscardCriterion& a, const DiscardCriterion& b) { return (a.quantity < b.quantity); } );

    EventQuantity lastqty = EventQuantity_None;

    // `perm` is a sort permutation and also maintains a deletion flag for each element
    perm.resize(nevents);  // `resize` and assign actually seems to be faster than `reserve` and `emplace_back`
    for (size_t i = 0; i < nevents; ++i)
      perm[i] = (i << 1);  // LSB is deletion flag; compute `(perm[...] >> 1)` to convert to an index

    for (auto& cut : discard_criteria_) {
      auto qty = cut.quantity;
      const auto& values = data_[qty];

      if (qty != lastqty) {
        std::sort( std::begin(perm), std::end(perm), [&values](size_t a, size_t b) { return values[a >> 1] < values[b >> 1]; } );  // always sorts ascending
        lastqty = qty;
      }

      size_t start;
      if (cut.is_percentile) {
        start = static_cast<size_t>(std::round(perm.size() * cut.value));

        cut.value = values[((start < perm.size()) ? perm[start] : perm.back()) >> 1];  // doesn't account for manual inclusions or first/last edges, but cuts shouldn't be that close to the ends anyway
        cut.is_percentile = false;
      }
      else {
        start = binary_search_ascending_for_discard(perm, values, cut.value);

        if (!(cut.greater_than) && start < perm.size() && values[perm[start] >> 1] == cut.value)
          ++start;  // exact matches are always discarded, for both `<` and `>`
      }

      if (cut.greater_than) {
        for (auto it = std::begin(perm) + static_cast<std::ptrdiff_t>(start); it != std::end(perm); ++it)
          *it |= 1;  // set deletion flag for this event
      }
      else {
        for (auto it = std::begin(perm), itend = it + static_cast<std::ptrdiff_t>(start); it != itend; ++it)
          *it |= 1;  // set deletion flag for this event
      }
    }

    std::sort(perm.begin(), perm.end());  // numerically sort `perm` to order deletion flags by the order in which events were recorded

    for (auto& field : data_) {
      remove_for_discard(perm, field.second);
      nevents = field.second.size();  // easy way to count how many events were kept; should be the same for all `field`s
    }

    perm.resize(0);  // empty `perm` so it can be reused (without reallocation) below
  }

  // resolve percentile bin edges

  for (auto& axis : axes_) {
    if (!axis.is_percentile)
      continue;

    axis.is_percentile = false;

    if (nevents == 0) {
      for (auto& edge : axis.edges)
        edge = 0.;
      continue;
    }

    // now `perm` is just a sort permutation (LSB is not special)
    if (perm.size() == 0) {
      perm.resize(nevents);  // `resize` and assign actually seems to be faster than `reserve` and `emplace_back`
      for (size_t i = 0; i < nevents; ++i)
        perm[i] = i;
    }

    const auto& values = data_[axis.quantity];

    std::sort( std::begin(perm), std::end(perm), [&values](size_t a, size_t b) { return values[a] < values[b]; } );

    for (auto& edge : axis.edges)
      edge = values[perm[std::min(static_cast<size_t>(std::round(nevents * edge)), nevents - 1)]];

    static constexpr double SMALL = 1E-6;

    if (axis.has_first_inclusion) {
      if (axis.ascending)
        axis.edges.front() = std::min(axis.edges.front(), axis.first_inclusion);
      else axis.edges.front() = std::max(axis.edges.front(), axis.first_inclusion);
    }
    else if (!axis.ascending) {
      if (std::abs(axis.edges.front()) < SMALL)
        axis.edges.front() = SMALL;
      else axis.edges.front() += (std::abs(axis.edges.front()) * SMALL);
    }

    if (axis.has_last_inclusion) {
      if (axis.ascending)
        axis.edges.back() = std::max(axis.edges.back(), axis.last_inclusion);
      else axis.edges.back() = std::min(axis.edges.back(), axis.last_inclusion);
    }
    else if (axis.ascending) {
      if (std::abs(axis.edges.back()) < SMALL)
        axis.edges.back() = SMALL;
      else axis.edges.back() += (std::abs(axis.edges.back()) * SMALL);
    }
  }

  // bin all pending events into the finalized bins

  for (size_t i = 0; i < nevents; ++i) {
    bin_event( [this, i](EventQuantity qty) -> double { return data_[qty][i]; } );
  }

  data_.clear();
  numevents_pending_ = 0;
  delay_binning_ = false;
}


std::pair<double, bool> EventBinnerHelper::parse_cutoff(const std::string& str) {
  double d;

  if (!boost::ends_with(str, "%")) {
    d = std::stod(str);
    if (std::isnan(d))
      throw std::invalid_argument{"quantity value must be a real number"};

    return std::make_pair(d, false);
  }

  if (boost::ends_with(str, "_%"))
    d = std::stod(str.substr(0, str.size() - 2));
  else d = 100. - std::stod(str.substr(0, str.size() - 1));

  if (d < 0. || d > 100. || !std::isfinite(d))
    throw std::invalid_argument{"percentile must be a real number from 0 to 100"};

  return std::make_pair(d / 100., true);
}

void EventBinnerHelper::parse_edges(std::string&& paramstr) {
  boost::char_separator<char> sep(",", "", boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tokens(paramstr, sep);

  auto it = std::begin(tokens);
  if (it == std::end(tokens))
    throw std::invalid_argument{"bin edge list must begin with a quantity label"};

  EventQuantity qty = EventQuantityList::lookup(*it);

  std::vector<double> edges;
  bool pct = false;

  double inclusion1 = 0., inclusion2 = 0.;
  bool has_inclusion1 = false, has_inclusion2 = false;

  while (++it != std::end(tokens)) {
    auto edge = parse_cutoff(*it);

    if (edge.second) {  // is a percentile
      if (!pct && edges.size() == 1) {
        inclusion1 = edges.back();
        has_inclusion1 = true;
        edges.pop_back();
      }
      else if ((!pct && edges.size() != 0) || has_inclusion2)
        throw std::invalid_argument{"bin edges cannot be a mixture of values and percentiles"};

      edges.push_back(edge.first);
      pct = true;
    }
    else {  // not a percentile
      if (pct) {
        if (!has_inclusion2) {
          inclusion2 = edge.first;
          has_inclusion2 = true;
        }
        else throw std::invalid_argument{"bin edges cannot be a mixture of values and percentiles"};
      }
      else edges.push_back(edge.first);
    }
  }

  if (!pct) {
    binner_.add_axis(qty, std::move(edges));
  }
  else if (has_inclusion1) {
    if (has_inclusion2)
      binner_.add_axis_percentile(qty, std::move(edges), inclusion1, inclusion2);
    else binner_.add_axis_percentile(qty, std::move(edges), inclusion1, false);
  }
  else if (has_inclusion2)
    binner_.add_axis_percentile(qty, std::move(edges), inclusion2, true);
  else binner_.add_axis_percentile(qty, std::move(edges));

  require_quantity_callback_(qty);
}

void EventBinnerHelper::parse_discard(std::string&& inequality) {
  boost::char_separator<char> sep("", "<>", boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tokens(inequality, sep);

  std::vector<std::string> tokenlist(std::begin(tokens), std::end(tokens));
  if (tokenlist.size() != 3 || (tokenlist[1] != "<" && tokenlist[1] != ">"))
    throw std::invalid_argument{"discard directive requires a strict inequality"};

  EventQuantity qty = EventQuantityList::lookup(tokenlist.front());
  bool gt = (tokenlist[1] == ">");
  auto&& cutoff = parse_cutoff(tokenlist.back());

  if (cutoff.second)
    binner_.add_discard_percentile(qty, gt, cutoff.first);
  else binner_.add_discard(qty, gt, cutoff.first);

  require_quantity_callback_(qty);
}

void EventBinnerHelper::parse_directives(const std::string& directives) {
  boost::char_separator<char> sep(":", "", boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tokens(directives, sep);

  for (auto it = ++std::begin(tokens); it != std::end(tokens); ++it) {
    const auto& token = *it;
    auto x = token.find('=');
    std::string what = token.substr(0, x);

    if (what == "edges") {
      if (x == std::string::npos)
        throw std::invalid_argument{"edges directive requires a list"};

      parse_edges(token.substr(x + 1));
    }
    else if (what == "subsample") {
      if (x == std::string::npos)
        throw std::invalid_argument{"subsample directive requires a nonnegative integer"};

      binner_.set_subsample(std::stoul(token.substr(x + 1)));
    }
    else if (what == "discard") {
      if (x == std::string::npos)
        throw std::invalid_argument{"discard directive requires an inequality"};

      parse_discard(token.substr(x + 1));
    }
    else {
      unhandled_directive_callback_(token);
    }
  }
}

void EventBinnerHelper::configure(EventBinner& binner, const std::string& directives,
  std::function<void(const std::string&)> unhandled_directive_callback,
  std::function<void(EventQuantity)> require_quantity_callback)
{
  EventBinnerHelper helper(binner, unhandled_directive_callback, require_quantity_callback);

  helper.parse_directives(directives);

  binner.start();
}

}  // namespace trento

