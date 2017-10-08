#include <iostream>


#ifdef STANDALONE
#include "Event.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#endif

Event::Event(uint32_t run, uint32_t event): run_(run), event_(event) {
}

void Event::fillHit(uint32_t id, double energy) {
    auto itr_insert = hits_.emplace(id, energy);
    // if id already inserted (noise and hit), add the energy
    if(!itr_insert.second) itr_insert.first->second += energy;
}

void Event::fillGenEn(uint32_t id, double energy) {
    generated_energy_.emplace(id, energy);
}

void Event::fillGenEta(uint32_t id, double eta) {
    generated_eta_.emplace(id, eta);
}

void Event::fillGenPhi(uint32_t id, double phi) {
    generated_phi_.emplace(id, phi);
}

void Event::fillPDGid(uint32_t id, int PDGid) {
    pdgid_.emplace(id, PDGid);
}

void Event::fillThick(uint32_t id, int thick) {
    if (!(thick_.find(id) != thick_.end()))
        thick_.emplace(id, thick);
}

void Event::fillCells(uint32_t id, Cell& cell) {
    // if the cell is already in map, don't add it
    if (!(cells_.find(id) != cells_.end()))
        cells_.emplace(id, cell);
}

void Event::setnPart(uint32_t part) {
    npart_ = part;
}

void Event::clear() {
    hits_.clear();
}


int Event::getLayerFromId(int cell_id) {

    int len = std::to_string(cell_id).length();

    return atoi(std::to_string(cell_id).substr(len -2, len).c_str());
}
