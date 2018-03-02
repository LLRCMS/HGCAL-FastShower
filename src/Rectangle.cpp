
#ifdef STANDALONE
#include "Rectangle.h"
#else
#include "HGCalSimulation/FastShower/interface/Rectangle.h"
#endif

//
// Contructors
//

Rectangle::
Rectangle(const Point& topLeft, const Point& bottomRight)
{
  this->topLeft_ = topLeft;
  this->bottomRight_ = bottomRight;
}

Rectangle::~Rectangle() { }

//
// Methods
//

bool Rectangle::
contains(const Point& p)
{
  return p.x >= topLeft_.x
      && p.x <= bottomRight_.x
      && p.y >= bottomRight_.y
      && p.y <= topLeft_.y;
}

const Point Rectangle::
getCenter()
{
  float dx = bottomRight_.x - topLeft_.x;
  float dy = topLeft_.y - bottomRight_.y;

  return Point(topLeft_.x + (dx / 2), bottomRight_.y + (dy / 2));
}

//
// Getters
//

const Point& Rectangle::
getTopLeft()
{
  return topLeft_;
}

const Point& Rectangle::
getBottomRight()
{
  return bottomRight_;
}
