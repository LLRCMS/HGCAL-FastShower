
#include "TMath.h"
#ifdef STANDALONE
#include "Cell.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#endif

#include <vector>


uint32_t Cell::id(uint16_t i, uint16_t j, uint16_t k) {
    // id format = iiijjjkk
    return i*100000 + j*100 + k;
}


Cell::Cell(TVectorD&& position,
           std::vector<TVectorD>&& vertices,
           int i_index,
           int j_index,
           int k_index):
            position_(std::move(position)),
            vertices_(std::move(vertices))
{

  i_index_ = (int16_t)i_index;
  j_index_ = (int16_t)j_index;
  k_index_ = (int16_t)k_index;
  id_ = Cell::id(i_index_, j_index_, k_index_);
  x_ = this->getPosition()(0);
  y_ = this->getPosition()(1);
  z_ = this->getPosition()(2);
}

double Cell::getX() const {
    return x_;
}

double Cell::getY() const{
    return y_;
}

double Cell::getZ() const{
    return z_;
}

int Cell::getIIndex() const{
    return i_index_;
}

int Cell::getJIndex() const{
    return j_index_;
}

int Cell::getLayer() const{
    return k_index_;
}

uint32_t Cell::getId() const{
    return id_;
}
