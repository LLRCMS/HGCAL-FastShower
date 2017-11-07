
#ifdef STANDALONE
#include "Point.h"
#else
#include "HGCalSimulation/FastShower/interface/Point.h"
#endif

Point::
Point(float x, float y)
{
    this->x = x;
    this->y = y;
}