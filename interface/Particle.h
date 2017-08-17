#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Parameters.h"

class Particle{

public:
    Particle(const Parameters&);
    ~Particle();

    double getEIncident();
    double getXIncident();
    double getYIncident();

private:
    double e_inc_;
    double x_incident;
    double y_incident;

    int nhits_;
    double denrj_;

    const Parameters& parameters_;

}

#endif