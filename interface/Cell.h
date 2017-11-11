
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
    Cell(TVectorD&&, std::vector<TVectorD>&&, double, int, int, int);

    ~Cell() {};

    // getters
    const TVectorD& getPosition() const {return position_;}
    const std::vector<TVectorD>& getVertices() const {return vertices_;}

    double getOrientation() const {return orientation_;}

    double getX() const;
    double getY() const;
    double getZ() const;
    int getIIndex() const;
    int getJIndex() const;
    int getLayer() const;
    uint32_t getId() const;

    bool isFullCell() const {return orientation_==90.;}
    bool isHalfCell() const {return int(orientation_)%60==0;}

  private:
    TVectorD position_; // centre position in absolute coordinates
    std::vector<TVectorD> vertices_; // vertices positions in absolute coordinates
    // hopefully 16 bits are enough
    uint16_t i_index_;
    uint16_t j_index_;
    uint16_t k_index_;
    // | j_index  | i_index | k_index
    uint32_t id_;

    double orientation_; // orientation for halh cells
};

#endif
