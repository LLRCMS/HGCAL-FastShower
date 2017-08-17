

#ifdef STANDALONE
#include "Event.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#endif


Event::Event(uint32_t run, uint32_t event): run_(run), event_(event) {
}

void Event::fillHit(uint32_t id, double energy) {
    auto itr_insert = hits_.emplace(id, energy);
    // if id already inserted, add the energy
    if(!itr_insert.second) itr_insert.first->second += energy;
}

void Event::fillGenEn(uint32_t id, double energy) {
    auto itr_insert = generated_energy_.emplace(id, energy);
    // if id already inserted, add the energy
    if(!itr_insert.second) itr_insert.first->second += energy;
}

void Event::fillGenEta(uint32_t id, double eta) {
    auto itr_insert = generated_eta_.emplace(id, eta);
    // if id already inserted, add the eta
    if(!itr_insert.second) itr_insert.first->second += eta;
}

void Event::fillGenPhi(uint32_t id, double phi) {
    auto itr_insert = generated_phi_.emplace(id, phi);
    // if id already inserted, add the phi
    if(!itr_insert.second) itr_insert.first->second += phi;
}

void Event::fillPDGid(uint32_t id, int PDGid) {
    auto itr_insert = pdgid_.emplace(id, PDGid);
    // if id already inserted, add the phi
    if(!itr_insert.second) itr_insert.first->second += PDGid;
}

void Event::fillThick(uint32_t id, int thick) {
    auto itr_insert = thick_.emplace(id, thick);
    // if id already inserted, add the phi
    if(!itr_insert.second) itr_insert.first->second += thick;
}

void Event::setnPart(uint32_t part) {
    npart_ = part;
}

void Event::clear() {
    hits_.clear();
}
