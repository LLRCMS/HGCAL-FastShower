#ifndef __HGCalSimulation_FastShower_Parameters_h__
#define __HGCalSimulation_FastShower_Parameters_h__

#include <map>
#include <string>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>

#include <vector>

const int NB_LAYERS = 52;
const int NB_LAYERS_EE = 28;
const int NB_SI_THICKNESS = 3;

class Parameters {
  public:
    struct General {
      General();
      unsigned events;
      std::vector<int> part_type;
      bool debug;
      std::string output_file;
    };

    struct Geometry {
      enum class Type {Hexagons, Triangles, External, Undefined};
      const static std::map<std::string, Type> type_map_;
      Geometry();
      Type type;
      int layer;
      std::vector<double> layers_z;
      // internal infinite geometries
      double small_cell_side;
      double large_cell_side;
      int EE_limit_layer;
      int FH_limit_layer;
      int BH_limit_layer;
      double limit_first_zone;
      double limit_second_zone;
      double eta_min;
      double eta_max;
      double phi_min;
      double phi_max;
      // external geometries
      std::string file;
    };

    struct Shower {
      Shower();
      double moliere_radius;
      std::map<std::string, double> transverse_parameters_electro;
      std::map<std::string, double> transverse_parameters_hadro;
      // std::vector<double> layers_energy;
      std::map<int, std::vector<double>> map_layers_energy;
      std::map<int, std::vector<double>> map_alpha;
    };

    struct Generation {
      enum class GenType {External, Personnal, Undefined};
      const static std::map<std::string, GenType> calib_type_map_;
      Generation();
      GenType gentype;
      bool fluctuation;
      bool fluctuation_energy;
      std::string gun_type;
      double energy;
      double E_range_min;
      double E_range_max;
      // own values for calibration
      double sampling;
      int number_of_hits_per_gev;
      bool noise;
      double noise_sigma;
      double incident_eta;
      double incident_phi;
      bool eta_fluctuation;
      bool phi_fluctuation;
      double eta_range_max;
      double eta_range_min;
      double phi_range_max;
      double phi_range_min;
      // calibration values read in file
      std::string calib_file;
    };

    struct Display {
      Display();
      unsigned events;
      int layer;
    };

public:
  Parameters();
  ~Parameters() {};

  void read(const std::string&);
  void print() const;

  const General& general() const {return general_;}
  const Geometry& geometry() const {return geometry_;}
  const Shower& shower() const {return shower_;}
  const Generation& generation() const {return generation_;}
  const Display& display() const {return display_;}

private:
  void fillGeneral(const boost::python::dict& dict);
  void fillGeometry(const boost::python::dict& dict);
  void fillShower(const boost::python::dict& dict);
  void fillGeneration(const boost::python::dict& dict);
  void fillDisplay(const boost::python::dict& dict);

  // Converters to std objects
  template<typename T>
  std::vector<T> toStdVector(const boost::python::list& plist) {
      boost::python::stl_input_iterator<T> begin(plist), end;
      return std::vector<T>(begin,end);
  }
  template<typename Key, typename Value>
  std::map<Key,Value> toStdMap(const boost::python::dict& pdict) {
    std::map<Key, Value> output;
    boost::python::list keys = pdict.keys();

    for(unsigned i=0; i<len(keys);++i) {
      Key key = boost::python::extract<Key>(keys[i]);
      Value value = boost::python::extract<Value>(pdict[key]);
      output.emplace(key, value);
    }
    return output;
  }

  template<typename Key, typename Value>
  std::map<Key, std::vector<Value>> toStdMapVector(const boost::python::dict& pdict) {

    std::map<Key, std::vector<Value>> output;
    boost::python::list keys = pdict.keys();

    for(unsigned i=0; i<len(keys);++i) {
      Key key = boost::python::extract<Key>(keys[i]);
      boost::python::list value_list = boost::python::extract<boost::python::list>(pdict[key]);
      std::vector<Value> vector_value = toStdVector<Value>(value_list);
      output.emplace(key, vector_value);
    }
    return output;
  }

  General general_;
  Geometry geometry_;
  Shower shower_;
  Generation generation_;
  Display display_;

};

#endif
