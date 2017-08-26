
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
    Cell(TVectorD&&, std::vector<TVectorD>&&, int, int, int);

    ~Cell() {};

    // getters
    const TVectorD& getPosition() const {return position_;}
    const std::vector<TVectorD>& getVertices() const {return vertices_;}


    double getX() const;
    double getY() const;
    double getZ() const;
    int getIIndex();
    int getJIndex();
    int getLayer() const;
    uint32_t getId() const;

  private:
    TVectorD position_; // centre position in absolute coordinates
    std::vector<TVectorD> vertices_; // vertices positions in absolute coordinates
    uint16_t i_index_;
    uint16_t j_index_;
    uint16_t k_index_;
    uint32_t id_;

    double x_;
    double y_;
    double z_;
};

#endif
