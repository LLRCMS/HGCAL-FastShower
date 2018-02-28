

#include "TMath.h"
#ifdef STANDALONE
#include "Cell.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#endif

#include <vector>


uint32_t Cell::id(uint16_t i, uint16_t j, uint16_t k) {
    // id format = iiijjjkk
    // return i*100000 + j*100 + k;
  uint32_t id = 0;
  id |= k & 0x3F;
  id |= ((j & 0x1FFF) << 6);
  id |= ((i & 0x1FFF) << 19);

return id;
}


Cell::Cell(TVectorD&& position, std::vector<TVectorD>&& vertices, double orientation, uint16_t i_index, uint16_t j_index, uint16_t k_index):
  position_(std::move(position)),
  vertices_(std::move(vertices)),
  orientation_(orientation)
{

  if(i_index > 0x1FFF || j_index > 0x1FFF ) {
    throw std::string("Cell index outside of 13bits integer limits");
  }
  i_index_ = (uint16_t)i_index;
  j_index_ = (uint16_t)j_index;
  k_index_ = (uint16_t)k_index;
  id_ = Cell::id(i_index_, j_index_, k_index_);
}

double Cell::getX() const {
    return this->getPosition()(0);
}

double Cell::getY() const{
    return this->getPosition()(1);
}

double Cell::getZ() const{
    return this->getPosition()(2);
}

uint16_t Cell::getIIndex() const{
    return i_index_;
}

uint16_t Cell::getJIndex() const{
    return j_index_;
}

uint16_t Cell::getLayer() const{
    return k_index_;
}

uint32_t Cell::getId() const{
    return id_;
}
