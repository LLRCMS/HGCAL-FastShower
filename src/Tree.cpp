
#include <iostream>

#ifdef STANDALONE
#include "Tree.h"
#else
#include "HGCalSimulation/FastShower/interface/Tree.h"
#endif

Tree::
Tree(const Rectangle& r, int levels)
{
  rectangle_ = r;
  // cells_ = nullptr;
  cells_.empty();

  subdivide(levels);
}

Tree::
~Tree()
{
  if(cells_.size() != 0) {
    cells_.empty();
  }

  nw_.reset();
  sw_.reset();
  ne_.reset();
  se_.reset();
}

//
// Methods
//

bool Tree::
empty()
{
  if (cells_.size() == 0) return true;
  else return false;
}

void Tree::
addCell(const Cell* c)
{
  cells_.push_back(c);
}

int Tree::
countCells()
{
  if(cells_.size() != 0) {
    // std::cout << "This leaf has " << cells_->size() << " elements" << std::endl;
    return cells_.size();
  }

  int cellCount = 0;

  cellCount += nw_->countCells();
  cellCount += sw_->countCells();
  cellCount += ne_->countCells();
  cellCount += se_->countCells();

  return cellCount;
}

const std::vector< const Cell*> Tree::
getCells()
{
  return cells_;
}

Tree* Tree::
getLeaf(float x, float y)
{
  const Point p(x, y);
  return getLeaf(p);
}

Tree* Tree::
getLeaf(const Point p)
{
  if(cells_.size() != 0) return this;

  if(nw_->rectangle_.contains(p)) {
    // std::cout << "Going to subtree nw_" << std::endl;
    return nw_->getLeaf(p);
  }

  if(sw_->rectangle_.contains(p)) {
    // std::cout << "Going to subtree sw" << std::endl;
    return sw_->getLeaf(p);
  }

  if(ne_->rectangle_.contains(p)) {
    // std::cout << "Going to subtree ne" << std::endl;
    return ne_->getLeaf(p);
  }

  if(se_->rectangle_.contains(p)) {
    // std::cout << "Going to subtree se" << std::endl;
    return se_->getLeaf(p);
  }

  // std::cout << "Looking for (" << p.x << ", " << p.y << ") "<< std::endl;
  // std::cout << "topLeft is (" << rectangle_.getTopLeft().x << ", "
  //     << rectangle_.getTopLeft().y << ")" << std::endl;
  // std::cout << "BottomRight is (" << rectangle_.getBottomRight().x << ", "
  //     << rectangle_.getBottomRight().y << ")" << std::endl;
  throw std::string("Coordinates not found.");
}

void Tree::
subdivide(int remainingLevels)
{
  const Point center = rectangle_.getCenter();

  const Rectangle r_nw(rectangle_.getTopLeft(), center);

  const Rectangle r_sw(
    Point(rectangle_.getTopLeft().x, center.y),
    Point(center.x, rectangle_.getBottomRight().y)
  );

  const Rectangle r_ne(
    Point(center.x, rectangle_.getTopLeft().y),
    Point(rectangle_.getBottomRight().x, center.y)
  );

  const Rectangle r_se(center, rectangle_.getBottomRight());

  nw_ = std::make_unique<Tree>(r_nw, remainingLevels - 1);
  sw_ = std::make_unique<Tree>(r_sw, remainingLevels - 1);
  ne_ = std::make_unique<Tree>(r_ne, remainingLevels - 1);
  se_ = std::make_unique<Tree>(r_se, remainingLevels - 1);
}