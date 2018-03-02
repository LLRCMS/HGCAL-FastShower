

#include "TMath.h"

#ifdef STANDALONE
#include "ShowerShapeHexagon.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShapeHexagon.h"
#endif



double ShowerShapeHexagon::firstNeighboors(int i, int j, int k) {

  // todo: add protections for the grid boarders
  int index;
  double sum=0.;
  sum += enrjMap_.at(Cell::id(i+1,j,k));
  sum += enrjMap_.at(Cell::id(i,j+1,k));
  sum += enrjMap_.at(Cell::id(i-1,j+1,k));
  sum += enrjMap_.at(Cell::id(i-1,j,k));
  sum += enrjMap_.at(Cell::id(i,j-1,k));
  sum += enrjMap_.at(Cell::id(i+1,j-1,k));
  return sum;
  
}
