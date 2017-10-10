
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


class Geometry {

  // a tesselation of the plane with polygonal cells
  // by default the plane is at z=0 but can be shifted to a z position that
  // corresponds to a layer position from TP geometry 

  public:

    Geometry(const Parameters::Geometry& params):parameters_(params){}
    ~Geometry();

    void constructFromParameters(bool, int, int);
    void constructFromJson(bool, int);


    double* dimensions(double side);
    double** hexagonoffset(double side);
    double** triangleoffset(double side);
    double* derivative(double side, Parameters::Geometry::Type itype);
    double** dxdyFirstZone(double* xs, double* ys);
    double** dxdySecondZone(double* xs, double* ys);
    int* ijWindows(int zone, double* xs, double* ys, double side, Parameters::Geometry::Type itype);
    double* XYrPhi(int i, int j, double side, Parameters::Geometry::Type itype, double* xs, double* ys, double zone);

    bool isInCell(const TVectorD& position, const Cell& cell) const; // test if a point is within a cell
    const Cell* closestCell(double x, double y) const;

    // getters
    std::unordered_map<uint32_t, Cell> getCells() const {return cells_;}
    int getLayer() const {return klayer_;}
    double getZlayer() const {return zlayer_;}

    const std::unique_ptr<TH2Poly>& cellHistogram() const {return cell_histogram_;}
    void draw(const Parameters::Display& params);
    void print();



  private:

    void setLayer(int klayer);
    std::string setHgcalPart(int klayer);
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