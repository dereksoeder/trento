// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>

#include "event.h"
#include "eventbin.h"
#include "hdf5_utils.h"

namespace trento {

namespace {

int compute_nevent_width(const VarMap& var_map) {
  // Determine the required width (padding) of the event number.  For example if
  // there are 10 events, the numbers are 0-9 and so no padding is necessary.
  // However given 11 events, the numbers are 00-10 with padded 00, 01, ...
  auto nevents = var_map["number-events"].as<int>();
  return static_cast<int>(std::ceil(std::log10(nevents)));
}

template<typename T = EventQuantity>
std::vector<T> without_duplicates(const std::vector<T>& list) {
  std::set<T> seen;
  std::vector<T> ret;

  for (const auto& elem : list) {
    if (seen.insert(elem).second)
      ret.push_back(elem);
  }

  return ret;
}

}  // unnamed namespace

// These output writer classes are used by the Output class.

class SummaryTextWriter : public OutputWriter {
 public:
  SummaryTextWriter(
    const VarMap& var_map,
    std::ostream& os,
    std::vector<EventQuantity>& required_quantities)
  : os_(os),
    width_(compute_nevent_width(var_map)),
    columns_(var_map["columns"].as<EventQuantityList>().values)
  {
    required_quantities.insert(std::end(required_quantities), std::begin(columns_), std::end(columns_));
  }

  virtual void process(const Event& event) {
    using std::fixed;
    using std::setprecision;
    using std::setw;
    using std::scientific;

    // Write a nicely-formatted line of event properties.
    // using a `stringstream` buffer seems to slightly improve performance
    buf_.str("");
    buf_.clear();
    buf_ << setprecision(10);

    for (auto it = std::begin(columns_); it != std::end(columns_); ++it) {
      auto qty = *it;

      // some special handling to ensure output is entirely backward-compatible
      // TODO: change all `setw(W)` to `' ' << setw(W-1)` to ensure columns never run together
      switch (qty) {
       case EventNumber:
        if (it != std::begin(columns_)) os_ << ' ';  // ensure `nevent` is padded if it isn't being displayed first
        buf_ << setw(width_) << event.nevent();
        break;
       case EventImpactParameter:
        buf_ << setw(15) << fixed << event.impact_param();
        break;
       case EventNumParticipants:
        buf_ << setw(5) << event.npart();
        break;
       case EventNumCollisions:
        if (event.ncoll() > 0)  // should always be positive if 'ncoll' column is present
          buf_ << setw(6) << event.ncoll();
        break;
       case EventMultiplicity:
        buf_ << setw(18) << scientific << event.multiplicity() << fixed;
        break;
       default:
        if (EventQuantityList::is_integer(qty))
          buf_ << ' ' << setw(4) << event.get<int>(qty);
        else buf_ << setw(14) << event.get<double>(qty);
        break;
      }
    }

    buf_ << '\n';
    os_ << buf_.str();
  }

 private:
  std::ostream& os_;
  int width_;
  const std::vector<EventQuantity> columns_;

  std::stringstream buf_;
};

class CSVSummaryTextWriter : public OutputWriter {
 public:
  CSVSummaryTextWriter(
    const VarMap& var_map,
    std::ostream& os,
    const std::string& directives,
    std::vector<EventQuantity>& required_quantities)
  : os_(os),
    columns_(var_map["columns"].as<EventQuantityList>().values)
  {
    boost::char_separator<char> sep(":", "", boost::keep_empty_tokens);
    boost::tokenizer<boost::char_separator<char>> tokens(directives, sep);

    for (auto it = ++std::begin(tokens); it != std::end(tokens); ++it) {  // skips initial "csv:"
      const auto& token = *it;

      if (token == "header")
        header_ = true;
      else throw std::invalid_argument{"unrecognized directive"};
    }

    required_quantities.insert(std::end(required_quantities), std::begin(columns_), std::end(columns_));
  }

  virtual void start() {
    if (header_) {
      buf_.str("");
      buf_.clear();

      bool newrow = true;

      for (auto qty : columns_) {
        if (newrow) newrow = false; else buf_ << ',';
        buf_ << '"' << EventQuantityList::lookup(qty) << '"';
      }

      buf_ << '\n';
      os_ << buf_.str();
    }
  }

  virtual void process(const Event& event) {
    buf_.str("");
    buf_.clear();
    buf_ << std::setprecision(10);

    bool newrow = true;

    for (auto qty : columns_) {
      if (newrow) newrow = false; else buf_ << ',';

      if (EventQuantityList::is_integer(qty))
        buf_ << event.get<int>(qty);
      else buf_ << event.get<double>(qty);
    }

    buf_ << '\n';
    os_ << buf_.str();
  }

 private:
  std::ostream& os_;
  const std::vector<EventQuantity> columns_;

  std::stringstream buf_;
  bool header_ = false;
} ;

class GridVisualizerTextWriter : public OutputWriter {
 public:
  GridVisualizerTextWriter(std::ostream& os) :
    os_(os)
  { }

  virtual void process(const Event& event) {
    const auto& grid = event.reduced_thickness_grid();

    double maxelem = TINY;
    for (const auto& row : grid)
      for (auto elem : row)
        maxelem = std::max(maxelem, elem);

    buf_.str("");
    buf_.clear();

    for (auto rowiter = grid.rbegin(); rowiter != grid.rend(); ++rowiter) {  // print highest-`y` row first (origin is bottom-left corner)
      for (auto elem : *rowiter) {
        if (elem < TINY)
          buf_ << ' ';
        else buf_ << scale_[static_cast<size_t>(scale_.size() * elem / maxelem)];
      }
      buf_ << '\n';
    }

    os_ << buf_.str();
  }

 private:
  std::ostream& os_;
  std::stringstream buf_;

  static constexpr double TINY = 1e-12;

  const std::string scale_ { ".,:/0123456789*" };
};

class FullTextWriter : public OutputWriter {
 public:
  FullTextWriter(
    const VarMap& var_map,
    const fs::path& output_dir,
    std::vector<EventQuantity>& required_quantities)
  : output_dir_(output_dir),
    width_(compute_nevent_width(var_map)),
    header_(!var_map["no-header"].as<bool>()),
    columns_(var_map["columns"].as<EventQuantityList>().values)
  {
    // Text files are all written into a single directory.  Require the
    // directory to be initially empty to avoid accidental overwriting and/or
    // spewing files into an already-used location.  If the directory does not
    // exist, create it.
    if (fs::exists(output_dir_)) {
      if (!fs::is_empty(output_dir_)) {
        throw std::runtime_error{"output directory '" + output_dir_.string() +
                                 "' must be empty"};
      }
    } else {
      fs::create_directories(output_dir_);
    }

    required_quantities.insert(std::end(required_quantities), std::begin(columns_), std::end(columns_));
  }

  virtual void process(const Event& event) {
    // Open a numbered file in the output directory.
    // Pad the filename with zeros.
    std::ostringstream padded_fname{};
    padded_fname << std::setw(width_) << std::setfill('0') << event.nevent() << ".dat";
    fs::ofstream ofs{output_dir_ / padded_fname.str()};

    if (header_) {
      // Write a commented header of event properties as key = value pairs.
      ofs << std::setprecision(10);

      for (auto qty : columns_) {
        switch (qty) {
         case EventNumber:
          ofs << "# event "   << event.nevent() << '\n';
          break;
         case EventNumCollisions:
          // Output ncoll if calculated
          if (event.ncoll() > 0)
            ofs << "# ncoll = " << event.ncoll() << '\n';
          break;
         default:
          ofs << "# " << std::left << std::setw(5) << EventQuantityList::lookup(qty) << " = ";
          if (EventQuantityList::is_integer(qty))
            ofs << event.get<int>(qty);
          else ofs << event.get<double>(qty);
          ofs << '\n';
          break;
        }
      }
    }

    // Write IC profile as a block grid.  Use C++ default float format (not
    // fixed-width) so that trailing zeros are omitted.  This significantly
    // increases output speed and saves disk space since many grid elements are
    // zero.
    for (const auto& row : event.reduced_thickness_grid()) {
      auto&& iter = row.begin();
      // Write all row elements except the last with a space delimiter afterwards.
      do {
        ofs << *iter << ' ';
      } while (++iter != --row.end());

      // Write the last element and a linebreak.
      ofs << *iter << '\n';
    }
  }

 private:
  const fs::path& output_dir_;
  int width_;
  bool header_;
  const std::vector<EventQuantity> columns_;
};

#ifdef TRENTO_HDF5

/// Simple class to write many events to an HDF5 file.
class HDF5Writer : public OutputWriter {
 public:
  /// Prepare an HDF5 file for writing.
  HDF5Writer(
    const VarMap& var_map,
    const fs::path& filename,
    std::vector<EventQuantity>& required_quantities)
  : file_(validate(filename).string(), H5F_ACC_TRUNC),
    columns_(without_duplicates(var_map["columns"].as<EventQuantityList>().values))
  {
    required_quantities.push_back(EventNumber);
    required_quantities.insert(std::end(required_quantities), std::begin(columns_), std::end(columns_));
  }

  /// Write an event.
  virtual void process(const Event& event) {
    // Prepare arguments for new HDF5 dataset.

    // The dataset name is a prefix plus the event number.
    const std::string name{"event_" + std::to_string(event.nevent())};

    // Cache a reference to the event grid -- will need it several times.
    const auto& grid = event.reduced_thickness_grid();

    // Define HDF5 datatype and dataspace to match the grid.
    const auto& datatype = hdf5::type<Event::Grid::element>();
    std::array<hsize_t, Event::Grid::dimensionality> shape;
    std::copy(grid.shape(), grid.shape() + shape.size(), shape.begin());
    auto dataspace = hdf5::make_dataspace(shape);

    // Set dataset storage properties.
    H5::DSetCreatPropList proplist{};
    // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
    // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
    // sense to chunk this way, since chunks must be read contiguously and there's
    // no reason to read a partial grid.
    proplist.setChunk(shape.size(), shape.data());
    // Set gzip compression level.  4 is the default in h5py.
    proplist.setDeflate(4);

    // Create the new dataset and write the grid.
    auto dataset = file_.createDataSet(name, datatype, dataspace, proplist);
    dataset.write(grid.data(), datatype);

    // Write event attributes.
    for (auto qty : columns_) {
      switch (qty) {
       case EventNumber:
        // already recorded event number in dataset's name
        break;
       case EventNumParticipants:
        hdf5_add_scalar_attr(dataset, EventQuantityList::lookup(qty), event.npart());
        break;
       case EventNumCollisions:
        if (event.ncoll() > 0)
          hdf5_add_scalar_attr(dataset, EventQuantityList::lookup(qty), event.ncoll());
        break;
      default:
        if (EventQuantityList::is_integer(qty))
          hdf5_add_scalar_attr(dataset, EventQuantityList::lookup(qty), event.get<int>(qty));
        else hdf5_add_scalar_attr(dataset, EventQuantityList::lookup(qty), event.get<double>(qty));
        break;
      }
    }
  }

 private:
  /// Add a simple scalar attribute to an HDF5 dataset.
  template <typename T>
  void hdf5_add_scalar_attr(
      const H5::DataSet& dataset, const std::string& name, const T& value) {
    const auto& datatype = hdf5::type<T>();
    auto attr = dataset.createAttribute(name, datatype, H5::DataSpace{});
    attr.write(datatype, &value);
  }

  // lets us check the status of the file before calling `file_`'s constructor
  const fs::path& validate(const fs::path& filename) {
    if (fs::exists(filename) && !fs::is_empty(filename))
      throw std::runtime_error{"file '" + filename.string() +
                               "' exists, will not overwrite"};
    return filename;
  }

  /// Internal storage of the file object.
  H5::H5File file_;

  const std::vector<EventQuantity> columns_;
};

#endif  // TRENTO_HDF5

class BinnedSummaryTextWriter : public OutputWriter {
 public:
  BinnedSummaryTextWriter(
    std::ostream& os,
    const VarMap& var_map,
    const std::string& directives,
    std::vector<EventQuantity>& required_quantities)
  : os_(os),
    columns_(var_map["columns"].as<EventQuantityList>().values),
    binner_(columns_)
  {
    required_quantities.insert(std::end(required_quantities), std::begin(columns_), std::end(columns_));

    EventBinnerHelper::configure(
      binner_,
      directives,
      [this](const std::string& token) {
        if (token == "noextent")
          show_extents_ = false;
        else if (token == "nocount")
          show_counts_ = false;
        else if (token == "nostdev")
          show_stdevs_ = false;
        else throw std::invalid_argument{"unrecognized directive"};
      },
      [&required_quantities](EventQuantity qty) {
        required_quantities.push_back(qty);
      }
    );
  }

  virtual void process(const Event& event) {
    binner_.process(event);
  }

  virtual void finish() {
    binner_.finish();

    os_ << std::setprecision(10);

    for (const auto& binresult : binner_) {
      bool newrow = true;

      if (show_extents_) {
        for (const auto& extent : binresult.extents) {
          if (newrow) newrow = false; else os_ << ',';
          os_ << extent.first << ',' << extent.second;
        }
      }

      if (show_counts_) {
        if (newrow) newrow = false; else os_ << ',';
        os_ << binresult.count;
      }

      for (const auto& result : binresult.data) {
        if (newrow) newrow = false; else os_ << ',';
        os_ << result.first;
        if (show_stdevs_) os_ << ',' << result.second;
      }

      os_ << '\n';
    }
  }

 private:
  std::ostream& os_;
  const std::vector<EventQuantity> columns_;
  EventBinner binner_;

  bool show_extents_ = true;
  bool show_counts_ = true;
  bool show_stdevs_ = true;
};

Output::Output(const VarMap& var_map) {
  bool quiet = var_map["quiet"].as<bool>();

  std::vector<std::string> outputdirectives;
  if (var_map.count("event-writer"))
    outputdirectives = var_map["event-writer"].as<std::vector<std::string>>();
  else outputdirectives.push_back("");  // default output writer only

  for (const auto& directive : outputdirectives) {
    auto&& writername = directive.substr(0, directive.find(':'));

    if (writername == "") {
      // Write to stdout unless the quiet option was specified.
      if (!quiet)
        writers_.emplace_back(
          std::unique_ptr<OutputWriter> { new SummaryTextWriter(var_map, std::cout, required_quantities_) }
        );
    }
    else if (writername == "csv" || writername == "CSV") {
      if (!quiet)
        writers_.emplace_back(
          std::unique_ptr<OutputWriter> { new CSVSummaryTextWriter(var_map, std::cout, directive, required_quantities_) }
        );
    }
    else if (writername == "gridvisualizer") {
      if (!quiet)
        writers_.emplace_back(
          std::unique_ptr<OutputWriter> { new GridVisualizerTextWriter(std::cout) }
        );
    }
    else if (writername == "binned") {
      if (!quiet)
        writers_.emplace_back(
          std::unique_ptr<OutputWriter> { new BinnedSummaryTextWriter(std::cout, var_map, directive, required_quantities_) }
        );
    }
    else throw std::invalid_argument{"unrecognized event writer"};
  }

  // Possibly write to text or HDF5 files.
  if (var_map.count("output")) {
    const auto& output_path = var_map["output"].as<fs::path>();
    if (hdf5::filename_is_hdf5(output_path)) {
#ifdef TRENTO_HDF5
      writers_.emplace_back(
        std::unique_ptr<OutputWriter> { new HDF5Writer(var_map, output_path, required_quantities_) }
      );
#else
      throw std::runtime_error{"HDF5 output was not compiled"};
#endif  // TRENTO_HDF5
    } else {
      writers_.emplace_back(
        std::unique_ptr<OutputWriter> { new FullTextWriter(var_map, output_path, required_quantities_) }
      );
    }
  }
}

void Output::start() const {
  for (const auto& writer : writers_)
    writer->start();
}

void Output::operator()(const Event& event) const {
  for (const auto& writer : writers_)
    writer->process(event);
}

void Output::operator()(int nevent, double impact_param, int ncoll, Event& event) const {
  // these external values are now stored in the event by Collider to be passed to the output writers;
  // this version of the function is retained for code compatibility but should not be used
  event.set_nevent(nevent);
  event.set_impact_parameter(impact_param);
  event.set_ncoll(ncoll);
  (*this)(event);
}

void Output::finish() const {
  for (const auto& writer : writers_)
    writer->finish();
}

}  // namespace trento
