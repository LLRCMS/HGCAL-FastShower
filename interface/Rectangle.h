#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include "Point.h"

class Rectangle {

private:
    Point* topLeft;
    Point* bottomRight;

public:
    Rectangle(Point*, Point*);
    ~Rectangle();

    bool contains(Point*);
    Point* getCenter();

    Point* getTopLeft();
    Point* getBottomRight();
};


#endif