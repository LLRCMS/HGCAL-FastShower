
#ifndef __HGCalSimulation_FastShower_Geometry_h__
#define __HGCalSimulation_FastShower_Geometry_h__

#include <string>
#include <vector>
#include <unordered_map>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH2Poly.h"
#ifdef STANDALONE
#include "Cell.h"
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif


const int NB_VTX_DIM = 4;
const int NB_OFFSET_MAX = 6;
const int NB_COORD_DIFF_MAX = 4; //triangle : up, up, down, down
const int NB_DERIV = 3;
const int NB_SEGMENT = 5;
const int NB_COORD_CART = 2;  //x and y
const int NB_PTS_WINDOWS = 4;  //imin, imax, jmin, jmax
const int NB_TOT_PTS_GEOMETRY = 9;
const int NB_COORD_CART_CYL = 4; //x, y, r, phi


class Geometry {

  // a tesselation of the plane with polygonal cells
  // by default the plane is at z=0 but can be shifted to a z position that
  // corresponds to a layer position from TP geometry 

  public:

    Geometry(const Parameters::Geometry& params):parameters_(params){};
    ~Geometry();

    void constructFromParameters(bool, int, int);
    void constructFromJson(bool, int);

    std::array<double, NB_VTX_DIM> dimensions(double side);
    std::array< std::array<double, NB_OFFSET_MAX>, NB_COORD_DIFF_MAX> hexagonoffset(double);
    std::array< std::array<double, NB_OFFSET_MAX>, NB_COORD_DIFF_MAX> triangleoffset(double);
    std::array<double, NB_DERIV> derivative(double side, Parameters::Geometry::Type);
    std::array< std::array<double, NB_SEGMENT>, NB_COORD_CART> dxdyFirstZone(std::array<double, NB_TOT_PTS_GEOMETRY>, std::array<double, NB_TOT_PTS_GEOMETRY>);
    std::array< std::array<double, NB_SEGMENT>, NB_COORD_CART>  dxdySecondZone(std::array<double, NB_TOT_PTS_GEOMETRY>, std::array<double, NB_TOT_PTS_GEOMETRY>);
    std::array< int, NB_PTS_WINDOWS> ijWindows(int, std::array<double, NB_TOT_PTS_GEOMETRY>, std::array<double, NB_TOT_PTS_GEOMETRY>, double, Parameters::Geometry::Type);
    std::array< double, NB_COORD_CART_CYL> XYrPhi(int, int, double, Parameters::Geometry::Type, std::array<double, NB_TOT_PTS_GEOMETRY>, std::array<double, NB_TOT_PTS_GEOMETRY>, double);

    bool isInCell(const TVectorD& position, const Cell& cell) const; // test if a point is within a cell
    const Cell* closestCell(double x, double y) const;

    // getters
    const std::unordered_map<uint32_t, Cell>& getCells() const {return cells_;}
    int getLayer() const {return klayer_;}
    double getZlayer() const {return zlayer_;}

    const std::unique_ptr<TH2Poly>& cellHistogram() const {return cell_histogram_;}
    void draw(const Parameters::Display& params);
    void print();



  private:

    void setLayer(int klayer);
    void setZlayer(double zlayer) {zlayer_ = zlayer;}

    std::unordered_map<uint32_t, Cell> cells_;
    std::unique_ptr<TH2Poly> cell_histogram_;
    Parameters::Geometry::Type itype_; // cell type
    int i_cell_first;
    int i_cell_second;
    int i_cell_third;
    int j_cell_first;
    int j_cell_second;
    int j_cell_third;
    int klayer_;
    const Parameters::Geometry& parameters_;
    double zlayer_;
};

#endif