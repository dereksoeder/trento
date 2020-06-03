// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015-2020 Jonah E. Bernhard, J. Scott Moreland, Derek Soeder
// MIT License

#include <stdexcept>
#include <boost/tokenizer.hpp>

#include "eventqty.h"

namespace trento {

/// Default quantities displayed in Trento output.
const EventQuantityList DefaultEventQuantityList {
  EventNumber,                  // nevent
  EventImpactParameter,         // b
  EventNumParticipants,         // npart
  EventMultiplicity,            // mult
  EventQuantity_Epsilon_n(2),   // e2
  EventQuantity_Epsilon_n(3),   // e3
  EventQuantity_Epsilon_n(4),   // e4
  EventQuantity_Epsilon_n(5)    // e5
} ;

/// Maps numerical identifiers to textual names for event quantities.
std::map<EventQuantity, std::vector<std::string>> EventQuantityLabels(
  []() {
    // note: the first label in each list is the one that gets returned by `lookup` and emitted by output writers
    std::map<EventQuantity, std::vector<std::string>> init {
      { EventNumber, { "nevent" } },
      { EventImpactParameter, { "b", "impact_param", "impact_parameter" } },
      { EventNumParticipants, { "npart" } },
      { EventNumCollisions, { "ncoll" } },
      { EventMultiplicity, { "mult" } },
      { EventFailuresA, { "failures_a" } },
      { EventFailuresB, { "failures_b" } },
      { EventGridCMx, { "cmx" } },
      { EventGridCMy, { "cmy" } }
    };

    for (int m = MinEccentricityHarmonic; m <= MaxEccentricityHarmonic; ++m)
      for (int n = MinEccentricityHarmonic; n <= MaxEccentricityHarmonic; ++n) {
        auto id = EventQuantity_Epsilon_mn(m,n);
        if (!init.count(id)) {
          auto& labels = init[id];  // `[]` operator creates entry with an empty vector, which `labels` will then reference

          if (m == n) labels.push_back("e" + std::to_string(n));
          labels.push_back("em" + std::to_string(m) + "n" + std::to_string(n));  // em2n2, em2n3, ... em8n7, em8n8
        }
      }

    return init;
  }()
);

EventQuantity EventQuantityList::lookup(const std::string& label) {
  // this could be faster, but it isn't time-critical and so doesn't justify the code complexity
  for (const auto& item : EventQuantityLabels)
    if (std::find(item.second.cbegin(), item.second.cend(), label) != item.second.cend())
      return item.first;

  throw std::invalid_argument{"unrecognized quantity"};
}

const std::string& EventQuantityList::lookup(EventQuantity id) {
  if (EventQuantityLabels.count(id) > 0) {
    const auto& labels = EventQuantityLabels[id];
    if (!labels.empty())
      return labels.front();
  }

  throw std::invalid_argument{"unrecognized quantity"};
}

void EventQuantityList::append(const std::string& token) {
  boost::char_separator<char> sep(",", "");
  boost::tokenizer<boost::char_separator<char>> subtokens(token, sep);

  for (const auto& label : subtokens)
    values.push_back(lookup(label));
}

/// Custom `boost::program_options` validator to pre-parse string argments into `EventQuantityList`s.
void validate(boost::any& v, const std::vector<std::string>& tokens, EventQuantityList *, int) {
  // see: https://stackoverflow.com/questions/3065109/can-boost-program-options-separate-comma-separated-argument-values
  // and: https://www.boost.org/doc/libs/1_73_0/doc/html/program_options/howto.html#id-1.3.31.6.7

  if (v.empty()) {
    v = boost::any(EventQuantityList());
  }

  auto p = boost::any_cast<EventQuantityList>(&v);

  for (const auto& token : tokens) {
    p->append(token);
  }
}

}  // namespace trento
