
#ifndef __HGCalSimulation_FastShower_Cell_h__
#define __HGCalSimulation_FastShower_Cell_h__

#include "TVectorD.h"

class Cell {

  public:
    static uint32_t id(uint16_t i, uint16_t j, uint16_t k);

  private:
    Cell() {}; // don't use default constructor


  public:
    // constructor for parametrized and full geometries
    // move position and vertices (the Cell object is the owner)
    // Cell(int, )
    Cell(TVectorD&&, std::vector<TVectorD>&&, double, uint16_t, uint16_t, uint16_t);

    ~Cell() {};

    // getters
    const TVectorD& getPosition() const {return position_;}
    const std::vector<TVectorD>& getVertices() const {return vertices_;}

    double getOrientation() const {return orientation_;}

    double getX() const;
    double getY() const;
    double getZ() const;
    uint16_t getIIndex() const;
    uint16_t getJIndex() const;
    uint16_t getLayer() const;
    uint32_t getId() const;

    bool isFullCell() const {return orientation_==90.;}
    bool isHalfCell() const {return int(orientation_)%60==0;}

  private:
    TVectorD position_; // centre position in absolute coordinates
    std::vector<TVectorD> vertices_; // vertices positions in absolute coordinates
    // on 32 bits, 6 are for the layer index, 26 are remaining for i and j
    uint16_t i_index_;
    uint16_t j_index_;
    uint16_t k_index_;
    // |  k_index  |  j_index  |  i_index  |
    // |32 ----- 19|18 ----- 6| 5 ----- 0|
    uint32_t id_;

    double orientation_; // orientation for halh cells
};

#endif
