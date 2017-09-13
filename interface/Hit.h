#ifndef HIT_H_
#define HIT_H_

#include "Parameters.h"
#include "ShowerParametrization.h"
#include "TVector.h"
#include "TRandom3.h"
#include "TMath.h"



class Hit {


public:

    Hit(const Parameters&);
    // Hit(int, int, double, double, double);
    ~Hit();

    double* hitPosition(int, int, double, double, double);
    TVector* getHitPosition();
    double getHitEnergy();

private:
    TRandom3 gun_;
    ShowerParametrization shower_;
    const Parameters& parameters_;


    TVectorD* position_;
    double* energy_;
};

#endif