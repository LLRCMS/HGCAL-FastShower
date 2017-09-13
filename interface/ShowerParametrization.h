#include <iostream>


#ifndef __HGCalSimulation_FastShower_ShowerParametrization_h__
#define __HGCalSimulation_FastShower_ShowerParametrization_h__

#ifdef STANDALONE
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif

using namespace std;

class ShowerParametrization {

    // Electromagnetic shower parametrization 
    // values from G. Grindhammer et al., hep-ex/0001020
    // adapted for HGCAL TP geometry

    public:

    ShowerParametrization() {}
    ShowerParametrization(const Parameters::Shower& params):
        moliereRadius_(params.moliere_radius),
        interaction_length_(params.interaction_length) {

        typedef map<int, std::vector<double>>::const_iterator MapIterator;
        for (MapIterator iter = params.map_layers_energy.begin(); iter != params.map_layers_energy.end(); iter++) {
            if(iter->second.size()!=52)
                throw std::string("The size of shower_layers_energy should be equals to the layer number");

            double total_weight=0.;
            typedef std::vector<double>::const_iterator ListIterator;
            for (ListIterator list_iter = iter->second.begin(); list_iter != iter->second.end(); list_iter++) {
                total_weight += *list_iter;
            }
            std::vector<double> new_list = iter->second;
            std::transform(new_list.begin(), new_list.end(), new_list.begin(),
                   std::bind2nd(std::divides<double>(),total_weight));
            layerProfile_[iter->first] = new_list;
        }

        for (MapIterator iter = params.map_alpha.begin(); iter != params.map_alpha.end(); iter++) {
            alpha_[iter->first] = iter->second;
        }

        // transverse parameters
        // exponential parameter set from TP studies, 90% containment at Rm for electro.
        // 95% containment at lambdaI for hadro
        r0layer_Electro_ = moliereRadius_/std::log(10.);
        r0layer_Hadro_ = interaction_length_/std::log(5.);
        const auto& transverse_params_electro = params.transverse_parameters_electro;
        // FIXME: find better way to implement transverse parameters
        a0_electro = transverse_params_electro.at("a0");
        a1_electro = transverse_params_electro.at("a1");
        a2_electro = transverse_params_electro.at("a2");

        const auto& transverse_params_hadro = params.transverse_parameters_hadro;
        a0_hadro = transverse_params_hadro.at("a0");
        a1_hadro = transverse_params_hadro.at("a1");
        a2_hadro = transverse_params_hadro.at("a2");

        tot_layer_Electro = 28.;
        tot_layer_Hadro = 52.;
    }

    ~ShowerParametrization() {}

    // average medium
    double getMoliereRadius() const {return moliereRadius_;}
    double getInteractionLength() const {return interaction_length_;}


    std::array<double,52>& getLayerProfile(int pdgid) {
        map<int, std::vector<double>>::const_iterator pos = layerProfile_.find(pdgid);
        if (pos == layerProfile_.end()) {
            cout << "pdgid doesn't exist"<<endl;
        } else {
            int index = 0;
            typedef std::vector<double>::const_iterator ListIterator;
            for (ListIterator list_iter = pos->second.begin(); list_iter != pos->second.end(); list_iter++){
                partLayerProfile_[index] = *list_iter;
                index++;
            }
        }
        return partLayerProfile_;
    }

    // transversal
    double r0_electro_(int klayer, int pdgid) {return (a0_electro + a1_electro*klayer + a2_electro*klayer*klayer)*r0layer_Electro_/tot_layer_Electro;}
    double r0_hadro_(int klayer, int pdgid) { return (a0_hadro + a1_hadro*klayer + a2_hadro*klayer*klayer)/tot_layer_Hadro;}

    // fluctuations
    std::array<double,3> spotEnergy(int pdgid) {
        map<int, std::vector<double>>::const_iterator pos = alpha_.find(pdgid);
        if (pos == alpha_.end()) {
            cout << "pdgid doesn't exist"<<endl;
        } else {
            int index = 0;
            typedef std::vector<double>::const_iterator ListIterator;
            for (ListIterator list_iter = pos->second.begin(); list_iter != pos->second.end(); list_iter++){
                alphasquare_[index] = (*list_iter)*(*list_iter);
                index++;
            }
        }
        return alphasquare_;}

    private:
        // average medium parameters
        std::array<double,52> radLength_;
        std::map<int, std::vector<double>> layerProfile_;
        std::array<double,52> partLayerProfile_;

        double moliereRadius_;   // in cm 
        double interaction_length_;  // in cm

        int tot_layer_Electro;
        int tot_layer_Hadro;

        // transverse parametrisation
        double r0layer_Electro_;
        double r0layer_Hadro_;
        double a0_hadro;
        double a1_hadro;
        double a2_hadro;
        double a0_electro;
        double a1_electro;
        double a2_electro;

        // for fluctuation
        std::map<int, std::vector<double>> alpha_;
        std::array<double,3> alphasquare_; // the stochastic coefficient in GeV^1/2
};

#endif
