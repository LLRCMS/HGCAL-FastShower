#include "Hit.h"

Hit::Hit(const Parameters& params): parameters_(params){
    position_ = nullptr;
    energy_ = nullptr;

    // hitPosition(layer, PDGid, e_inc, x_inc, y_inc);
}

Hit::~Hit(){
    if(position_ != nullptr) {
        delete position_;
    }
    if(energy_ != nullptr) {
        delete energy_;
    }
}

//
// Methods
//

// void Hit::hitCreation(int layer, int PDGid) {

//     // hit position on layer (EE, FH or BH). Energy attribution
//     int hit_pos = sqrt(x*x+y*y);
//     int thick = 0;
//     // double side;
//     if (hit_pos <= parameters_.geometry().limit_first_zone) {
//         real_energy = aShowerParametrization.spotEnergy(ip)[0];
//         thick = 100;
//         // side = parameters_.geometry().small_cell_side;
//     }
//     else if (hit_pos >= parameters_.geometry().limit_first_zone &&
//              hit_pos <= parameters_.geometry().limit_second_zone){
//         real_energy = aShowerParametrization.spotEnergy(ip)[1];
//         // side = parameters_.geometry().large_cell_side;
//         thick = 200;
//     }
//     else {
//         real_energy = denrj;
//         // side = parameters_.geometry().large_cell_side;
//         thick = 300;
//     }
//     energygen += real_energy;
// }




double * Hit::hitPosition(int layer, int PDGid, double e_inc, double x_inc, double y_inc) {

    // Shower parametrization
    ShowerParametrization aShowerParametrization(parameters_.shower());
    
    double r0_electro = aShowerParametrization.r0_electro_(layer, PDGid);
    double r0_hadro = aShowerParametrization.r0_hadro_(layer, PDGid);

    double f_em = std::log(e_inc)*0.1;

    double r_shower = 0;
    if (PDGid == 11 || PDGid == 22)
        r_shower = gun_.Exp(r0_electro); // exponential exp(-r/r0)
    else
        r_shower = f_em*gun_.Exp(r0_electro) + (1-f_em)*gun_.Exp(r0_hadro);

    double phi_shower = gun_.Rndm()*TMath::TwoPi();
    double x = r_shower*cos(phi_shower) + x_inc;
    double y = r_shower*sin(phi_shower) + y_inc;

    return new double[2] {x, y};
}