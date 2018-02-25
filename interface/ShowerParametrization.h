
#ifndef __HGCalSimulation_FastShower_ShowerParametrization_h__
#define __HGCalSimulation_FastShower_ShowerParametrization_h__

#ifdef STANDALONE
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif


class ShowerParametrization {

// Electromagnetic shower parametrization 
// values from G. Grindhammer et al., hep-ex/0001020
// adapted for HGCAL TP geometry

  public:

    ShowerParametrization() {}
    ShowerParametrization(const Parameters::Shower& params);

    ~ShowerParametrization() {}

    // average medium
    double getMoliereRadius() const {return moliereRadius_;}
    double getInteractionLength() const {return interaction_length_;}

    std::array<double, NB_LAYERS>& getLayerProfile(int pdgid);

    // transversal
    double r0_electro_(int klayer, int pdgid) {return (a0_electro + a1_electro*klayer + a2_electro*klayer*klayer)*r0layer_Electro_/tot_layer_Electro;}
    double r0_hadro_(int klayer, int pdgid) { return (a0_hadro + a1_hadro*klayer + a2_hadro*klayer*klayer)/tot_layer_Hadro;}

    // fluctuations
    double spotEnergy(int pdgid, int thickness);

private:
    // average medium parameters
    std::array<double, NB_LAYERS> radLength_;
    std::map<int, std::vector<double>> layerProfile_;
    std::array<double, NB_LAYERS> partLayerProfile_;

    double moliereRadius_;   // in cm 
    double interaction_length_;  // in cm

    int tot_layer_Electro;
    int tot_layer_Hadro;

    // transverse parametrisation
    double r0layer_Electro_;
    double a0_hadro;
    double a1_hadro;
    double a2_hadro;
    double a0_electro;
    double a1_electro;
    double a2_electro;

    // for fluctuation
    std::map<int, std::vector<double>> alpha_;
    double alphasquare_; // the stochastic coefficient in GeV^1/2
};

#endif
