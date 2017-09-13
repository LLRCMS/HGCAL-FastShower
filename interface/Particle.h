#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Parameters.h"

#include "TRandom3.h"

class Particle{

public:
    Particle(const Parameters&, int);
    ~Particle();

    double etaIncident(const Parameters&, TRandom3 );

    double getEIncident();
    double getXIncident();
    double getYIncident();

private:
    double e_inc_;
    double denrj_;

    TRandom3 gun_;

    int nhits_;

    const Parameters& parameters_;
    int particle_;

    double x_incident;
    double y_incident;

}

#endif