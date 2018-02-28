#ifndef __HGCalSimulation_FastShower_Rectangle_h__
#define __HGCalSimulation_FastShower_Rectangle_h__

#ifdef STANDALONE
#include "Point.h"
#else
#include "HGCalSimulation/FastShower/interface/Point.h"
#endif


class Rectangle {

  public:
    Rectangle(const Point& = 0, const Point& = 0);
    ~Rectangle();

    bool contains(const Point&);
    const Point getCenter();

    const Point& getTopLeft();
    const Point& getBottomRight();

  private:
    Point topLeft_;
    Point bottomRight_;
};

#endif