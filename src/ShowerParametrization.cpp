
#ifdef STANDALONE
#include "ShowerParametrization.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerParametrization.h"
#endif


ShowerParametrization::
ShowerParametrization(const Parameters::Shower& params): moliereRadius_(params.moliere_radius)
{
  for (auto& iter : params.map_layers_energy) {
    if(iter.second.size() != NB_LAYERS){
      throw std::string("The size of shower_layers_energy should be equals to the layer number");
    }

    double total_weight=0.;
    for (auto& list_iter : iter.second) total_weight += list_iter;

    std::vector<double> new_list = iter.second;
    std::transform(new_list.begin(), new_list.end(), new_list.begin(),
                   std::bind2nd(std::divides<double>(),total_weight));

    layerProfile_[iter.first] = new_list;
  }

  for (auto& iter : params.map_alpha) alpha_[iter.first] = iter.second;

  // transverse parameters
  // exponential parameter set from TP studies, 90% containment at Rm for electro.
  r0layer_Electro_ = moliereRadius_/std::log(10.);

  const auto& transverse_params_electro = params.transverse_parameters_electro;
  // FIXME: find better way to implement transverse parameters
  a0_electro = transverse_params_electro.at("a0");
  a1_electro = transverse_params_electro.at("a1");
  a2_electro = transverse_params_electro.at("a2");

  const auto& transverse_params_hadro = params.transverse_parameters_hadro;
  a0_hadro = transverse_params_hadro.at("a0");
  a1_hadro = transverse_params_hadro.at("a1");
  a2_hadro = transverse_params_hadro.at("a2");

  tot_layer_Electro = NB_LAYERS_EE;
  tot_layer_Hadro = NB_LAYERS;
}


std::array<double, NB_LAYERS>& ShowerParametrization::
getLayerProfile(int pdgid)
{

  map<int, std::vector<double>>::const_iterator pos = layerProfile_.find(pdgid);

  if (pos == layerProfile_.end()) cout << "pdgid doesn't exist"<<endl;
  else {
    int index = 0;
    for (auto& list_iter : pos->second){
      partLayerProfile_[index] = list_iter;
      index++;
    }
  }

  return partLayerProfile_;
}

std::array<double, NB_SI_THICKNESS> ShowerParametrization::
spotEnergy(int pdgid)
{

  map<int, std::vector<double>>::const_iterator pos = alpha_.find(pdgid);

  if (pos == alpha_.end()) {
      cout << "pdgid doesn't exist"<<endl;
  } else {
      int index = 0;
      for (auto& list_iter : pos->second){
          alphasquare_[index] = (list_iter)*(list_iter);
          index++;
      }
  }
  return alphasquare_;
}