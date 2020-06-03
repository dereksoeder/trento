// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef OUTPUT_H
#define OUTPUT_H

#include <memory>
#include <utility>
#include <vector>

#include "eventqty.h"
#include "fwd_decl.h"

namespace trento {

class OutputWriter {
 public:
  // needs to exist so derived writers will be destroyed when `Output::writers_` is destroyed
  virtual ~OutputWriter() = default;

  /// Writes any initial output before event processing begins.
  virtual void start() { }

  /// Handles the next event to be output.
  /// Implementations must not assume that 'event' will remain valid and unchanged after returning.
  virtual void process(const Event& event) = 0;

  /// Writes any final output after all events have been processed.
  /// May not be invoked in the event of an error.
  virtual void finish() { };
};

/// Simple interface for outputting event data.  Determines which output formats
/// to create based on the configuration and writes those formats when called.
class Output {
 public:
  /// Instantiate from the configuration.
  Output(const VarMap& var_map);

  /// Output any initial information before event processing begins.
  void start() const;

  /// Output event data.
  void operator()(const Event& event) const;

  /// Store external values in the event and output event data.
  void operator()(int nevent, double impact_param, int ncoll, Event& event) const;

  /// Output any final information after event processing ends.
  void finish() const;

  /// Get list of event quantities that output writers require to be computed.
  const std::vector<EventQuantity>& required_quantities() const
  { return required_quantities_; }

 private:
  /// Internal storage of output functions.
  std::vector<std::unique_ptr<OutputWriter>> writers_;

  /// Event quantities required by output writers.
  std::vector<EventQuantity> required_quantities_;
};

}  // namespace trento

#endif  // OUTPUT_H
