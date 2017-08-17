#ifndef HIT_H_
#define HIT_H_

#include "Parameters.h"
#include "ShowerParametrization.h"
#include "TVector.h"
#include "TRandom3.h"



class Hit {


public:

    Hit(int, int, double);
    ~Hit();

    TVector* hitPosition(int, int, double);
    TVector* getHitPosition();
    double getHitEnergy();

private:
    TRandom3 gun_;
    ShowerParametrization shower_;
    const Parameters& parameters_;


    TVector* position_;
    double* energy_;
};

#endif