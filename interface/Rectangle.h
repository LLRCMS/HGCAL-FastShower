#ifndef __HGCalSimulation_FastShower_Rectangle_h__
#define __HGCalSimulation_FastShower_Rectangle_h__

#ifdef STANDALONE
#include "Point.h"
#else
#include "HGCalSimulation/FastShower/interface/Point.h"
#endif


class Rectangle {

private:
  Point topLeft;
  Point bottomRight;

public:
  Rectangle(Point = 0, Point = 0);
  ~Rectangle();

  bool contains(Point);
  Point getCenter();

  Point getTopLeft();
  Point getBottomRight();
};

#endif