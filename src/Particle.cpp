#include "Particle.h"


Particle::Particle(int particle):particle_(particle) {
}

Particle::~Particle() {
}


// //
// Methods
// 

double Particle::etaIncident() {

    if (parameters_.generation().eta_fluctuation)
        return gun_.Uniform(parameters_.generation().eta_range_min,
                            parameters_.generation().eta_range_max);
    else
        return parameters_.generation().incident_eta;
}

double Particle::phiIncident() {

    if (parameters_.generation().phi_fluctuation)
        return gun_.Uniform(parameters_.generation().phi_range_min,
                            parameters_.generation().phi_range_max);
    else
        return parameters_.generation().incident_phi;
}

double energyIncident() {

    if (parameters_.generation().fluctuation_energy) {
        double rand_energy = gun_.Uniform(parameters_.generation().E_range_min,
                                          parameters_.generation().E_range_max);

        if (parameters_.generation().gun_type == "E")
            return rand_energy;
        else{
            return rand_energy*((exp(eta_incident*eta_incident)+1)/(exp(eta_incident*eta_incident)-1));
        }
    }
    else {
        if (parameters_.generation().gun_type == "E")
            return parameters_.generation().energy;
        else
            return parameters_.generation().energy*((exp(eta_incident*eta_incident)+1)/(exp(eta_incident*eta_incident)-1));
    }

}