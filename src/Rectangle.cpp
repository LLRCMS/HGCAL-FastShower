#include "Rectangle.h"
#include <iostream>

//
// Contructors
//

Rectangle::Rectangle(Point* topLeft, Point* bottomRight) {
    this->topLeft = topLeft;
    this->bottomRight = bottomRight;
}

Rectangle::~Rectangle() { }

//
// Methods
//

bool Rectangle::contains(Point* p) {
    return p->x >= topLeft->x
        && p->x <= bottomRight->x
        && p->y >= bottomRight->y
        && p->y <= topLeft->y;
}

Point* Rectangle::getCenter() {
    float dx = bottomRight->x - topLeft->x;
    float dy = topLeft->y - bottomRight->y;

    return new Point(topLeft->x + (dx / 2), bottomRight->y + (dy / 2));
}

//
// Getters
//

Point* Rectangle::getTopLeft() {
    return topLeft;
}

Point* Rectangle::getBottomRight() {
    return bottomRight;
}
