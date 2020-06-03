// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015-2020 Jonah E. Bernhard, J. Scott Moreland, Derek Soeder
// MIT License

#ifndef EVENTBIN_H
#define EVENTBIN_H

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "event.h"
#include "eventqty.h"
#include "fwd_decl.h"

namespace trento {

class EventBinnerHelper;

class EventBinner {
 public:
  EventBinner(const std::vector<EventQuantity>& columns);

  void set_subsample(size_t count);

  void add_axis(EventQuantity qty, std::vector<double>&& edges);

  void add_axis_percentile(EventQuantity qty, std::vector<double>&& edges);

  void add_axis_percentile(EventQuantity qty, std::vector<double>&& edges, double include_value, bool in_last_bin);

  void add_axis_percentile(EventQuantity qty, std::vector<double>&& edges, double include_first, double include_last);

  void add_discard(EventQuantity qty, bool greater_than, double value);

  void add_discard_percentile(EventQuantity qty, bool greater_than, double percentile);

  void start();

  void process(const Event& event);

  void finish();

  static void configure(EventBinner& binner, const std::string& directives);

  struct BinResult {
    std::vector<std::pair<double, double>> extents;  // {left edge, right edge} per axis
    std::vector<std::pair<double, double>> data;     // {mean, stdev} per column
    unsigned int count;                              // number of events that fell into this bin

    BinResult(size_t numcolumns) : extents(numcolumns) { }
  } ;

  // see: https://en.cppreference.com/w/cpp/iterator/iterator
  class bin_iterator : public std::iterator<std::input_iterator_tag, BinResult> {
   private:
    friend EventBinner;

    const EventBinner& binner_;
    size_t index_;

    explicit bin_iterator(const EventBinner& binner, size_t index = 0) : binner_(binner), index_(index) { }

   public:
    bin_iterator& operator++() { ++index_; return *this; }  // pre-increment
    bin_iterator operator++(int) { auto ret = *this; ++(*this); return ret; }  // post-increment
    bool operator==(bin_iterator other) const { return (index_ == other.index_); }
    bool operator!=(bin_iterator other) const { return !(*this == other); }

    BinResult operator*() const {
      BinResult ret(binner_.axes_.size());

      ret.count = binner_.bin_counts_[index_];

      if (ret.count > 0) {
        auto itsums = binner_.bin_sums_.cbegin() + static_cast<ptrdiff_t>(index_ * binner_.columns_.size());
        for (size_t n = binner_.columns_.size(); n != 0; ++itsums, --n) {
          const auto& binsums = *itsums;
          auto mean  = binsums.first / static_cast<double>(ret.count);
          auto stdev = ((ret.count > 1) ? std::sqrt((binsums.second / static_cast<double>(ret.count)) - (mean * mean)) : 0);
          ret.data.emplace_back( std::make_pair(mean, stdev) );
        }
      } else {
        ret.data.resize(binner_.columns_.size());  // use default values (zeroes)
      }

      auto tmpindex = index_;
      for (size_t i = ret.extents.size(); i != 0; ) {
        const auto& axis = binner_.axes_[--i];
        auto b = (tmpindex % (axis.edges.size() - 1));
        tmpindex /= (axis.edges.size() - 1);
        ret.extents.at(i) = std::make_pair(axis.edges[b], axis.edges[b + 1]);
      }

      return ret;
    }
  } ;

  bin_iterator begin() const {
    if (delay_binning_)
      throw std::logic_error{"cannot access bins until event percentiles have been resolved"};

    return bin_iterator(*this, 0);
  }

  bin_iterator end() const {
    return bin_iterator(*this, bin_counts_.size());
  }

 private:
  const std::vector<EventQuantity> columns_;
  std::map<EventQuantity, std::vector<double>> data_;

  struct BinAxis {
    BinAxis() = default;

    BinAxis(EventQuantity qty, std::vector<double>&& e, bool asc, bool pct) :
      quantity(qty), edges(std::move(e)), ascending(asc), is_percentile(pct),
      has_first_inclusion(false), has_last_inclusion(false)
    { }

    BinAxis(EventQuantity qty, std::vector<double>&& e, bool asc, double inclusion, bool is_last_bin) :
      quantity(qty), edges(std::move(e)), ascending(asc), is_percentile(true),
      has_first_inclusion(!is_last_bin), has_last_inclusion(is_last_bin),
      first_inclusion(is_last_bin ? double() : inclusion), last_inclusion(is_last_bin ? inclusion : double())
    { }

    BinAxis(EventQuantity qty, std::vector<double>&& e, bool asc, double inclusion1, double inclusion2) :
      quantity(qty), edges(std::move(e)), ascending(asc), is_percentile(true),
      has_first_inclusion(true), has_last_inclusion(true),
      first_inclusion(inclusion1), last_inclusion(inclusion2)
    { }

    EventQuantity quantity;
    std::vector<double> edges;
    bool ascending;
    bool is_percentile;
    bool has_first_inclusion;
    bool has_last_inclusion;
    double first_inclusion;
    double last_inclusion;
  } ;

  std::vector<BinAxis> axes_;

  struct DiscardCriterion {
    DiscardCriterion() = default;

    DiscardCriterion(EventQuantity qty, bool gt, double val, bool pct) :
      quantity(qty), greater_than(gt), value(val), is_percentile(pct)
    { }

    EventQuantity quantity;
    bool greater_than;
    double value;
    bool is_percentile;
  } ;

  std::vector<DiscardCriterion> discard_criteria_;

  size_t numevents_ = 0;
  size_t numevents_pending_ = 0;
  bool delay_binning_ = false;
  size_t subsample_ = 0;  // 0 = use all events to determine any bin percentiles

  bool should_discard(const Event& event) const;

  void resolve_percentiles();

  // should be multi-dimensional arrays, but the number of dimensions must be determined at runtime
  std::vector<unsigned int> bin_counts_;
  std::vector<std::pair<double, double>> bin_sums_;

  void bin_event(std::function<double(EventQuantity)> getter);
} ;

class EventBinnerHelper {
 public:
  static void configure(
    EventBinner& binner, const std::string& directives,
    std::function<void(const std::string&)> unhandled_directive_callback,
    std::function<void(EventQuantity)> require_quantity_callback);

  static void configure(EventBinner& binner, const std::string& directives) {
    configure(
      binner,
      directives,
      [](const std::string&) { throw std::invalid_argument{"unrecognized directive"}; },
      [](EventQuantity) { }
    );
  }

 private:
  EventBinnerHelper(
    EventBinner& binner,
    std::function<void(const std::string&)> unhandled_directive_callback,
    std::function<void(EventQuantity)> require_quantity_callback)
  : binner_(binner),
    unhandled_directive_callback_(unhandled_directive_callback),
    require_quantity_callback_(require_quantity_callback)
  { }

  static std::pair<double, bool> parse_cutoff(const std::string& str);

  void parse_edges(std::string&& paramstr);

  void parse_discard(std::string&& inequality);

  void parse_directives(const std::string& directives);

  EventBinner& binner_;
  const std::function<void(const std::string&)> unhandled_directive_callback_;
  const std::function<void(EventQuantity)> require_quantity_callback_;
} ;

}  // namespace trento

#endif  // EVENTBIN_H

