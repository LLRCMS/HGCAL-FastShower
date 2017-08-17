
#include "TMath.h"

#ifdef STANDALONE
#include "ShowerShapeTriangle.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShapeTriangle.h"
#endif



double ShowerShapeTriangle::firstNeighboors(int i, int j, int k) {

  // todo: add protections for the grid boarders
  int index;
  double sum=0.;
  if (imax_%2==0) { // upward triangles
    sum += enrjMap_.at(Cell::id(i+1,j, k));
    sum += enrjMap_.at(Cell::id(i-1,j, k));
    sum += enrjMap_.at(Cell::id(i+1,j-1, k));
  } else { // downward triangles
    sum += enrjMap_.at(Cell::id(i-1,j, k));
    sum += enrjMap_.at(Cell::id(i-1,j+1, k));
    sum += enrjMap_.at(Cell::id(i+1,j, k));
  }
  return sum;

}
