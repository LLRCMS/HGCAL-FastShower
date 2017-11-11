
#include <iostream>

#ifdef STANDALONE
#include "Tree.h"
#else
#include "HGCalSimulation/FastShower/interface/Tree.h"
#endif

Tree::
Tree(Rectangle r, int levels)
{
  rectangle = r;
  cells = nullptr;

  subdivide(levels);
}

Tree::
~Tree()
{
  if(cells != nullptr) {
    delete cells;
  }

  // delete rectangle;

  nw.reset();
  sw.reset();
  ne.reset();
  se.reset();
}

//
// Methods
//

bool Tree::
empty()
{
  if (cells == nullptr) return true;
  else return false;
}

void Tree::
addCell(Cell* c)
{
  cells->push_back(c);
}

int Tree::
countCells()
{
  if(cells != nullptr) {
    // std::cout << "This leaf has " << cells->size() << " elements" << std::endl;
    return cells->size();
  }

  int cellCount = 0;

  cellCount += nw->countCells();
  cellCount += sw->countCells();
  cellCount += ne->countCells();
  cellCount += se->countCells();

  return cellCount;
}

std::vector<Cell*>* Tree::
getCells()
{
  return cells;
}

Tree* Tree::
getLeaf(float x, float y)
{
  Point p = Point(x, y);
  return getLeaf(p);
}

Tree* Tree::
getLeaf(Point p)
{
  if(cells != nullptr) return this;

  if(nw->rectangle.contains(p)) {
    // std::cout << "Going to subtree nw" << std::endl;
    return nw->getLeaf(p);
  }

  if(sw->rectangle.contains(p)) {
    // std::cout << "Going to subtree sw" << std::endl;
    return sw->getLeaf(p);
  }

  if(ne->rectangle.contains(p)) {
    // std::cout << "Going to subtree ne" << std::endl;
    return ne->getLeaf(p);
  }

  if(se->rectangle.contains(p)) {
    // std::cout << "Going to subtree se" << std::endl;
    return se->getLeaf(p);
  }

  // std::cout << "Looking for (" << p.x << ", " << p.y << ") "<< std::endl;
  // std::cout << "topLeft is (" << rectangle.getTopLeft().x << ", "
  //     << rectangle.getTopLeft().y << ")" << std::endl;
  // std::cout << "BottomRight is (" << rectangle.getBottomRight().x << ", "
  //     << rectangle.getBottomRight().y << ")" << std::endl;
  throw std::string("Coordinates not found.");
}

void Tree::
subdivide(int remainingLevels)
{
  if(remainingLevels == 0) {
    cells = new std::vector<Cell*>;
    return;
  }

  Point center = rectangle.getCenter();

  Rectangle r_nw = Rectangle(rectangle.getTopLeft(), center);

  Rectangle r_sw = Rectangle(
    Point(rectangle.getTopLeft().x, center.y),
    Point(center.x, rectangle.getBottomRight().y)
  );

  Rectangle r_ne = Rectangle(
    Point(center.x, rectangle.getTopLeft().y),
    Point(rectangle.getBottomRight().x, center.y)
  );

  Rectangle r_se = Rectangle(center, rectangle.getBottomRight());

  nw = std::make_unique<Tree>(r_nw, remainingLevels - 1);
  sw = std::make_unique<Tree>(r_sw, remainingLevels - 1);
  ne = std::make_unique<Tree>(r_ne, remainingLevels - 1);
  se = std::make_unique<Tree>(r_se, remainingLevels - 1);
}