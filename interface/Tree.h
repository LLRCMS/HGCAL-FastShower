#ifndef __HGCalSimulation_FastShower_tree_h__
#define __HGCalSimulation_FastShower_tree_h__

#ifdef STANDALONE
#include "Cell.h"
#include "Rectangle.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#include "HGCalSimulation/FastShower/interface/Rectangle.h"
#endif

#include <vector>
# include <memory>

class Tree {

  public:
    Tree(Rectangle, int);
    ~Tree();

    bool empty();
    void addCell(const Cell*);
    int countCells();
    std::vector<const Cell*>* getCells();
    Tree* getLeaf(Point);
    Tree* getLeaf(float, float);

  private:
    //4 childrens
    std::unique_ptr<Tree> nw;
    std::unique_ptr<Tree> ne;
    std::unique_ptr<Tree> sw;
    std::unique_ptr<Tree> se;

  void subdivide(int);

  // Cells selected in the tree
  std::vector<const Cell*>* cells;

  protected:
  Rectangle rectangle;
};

#endif
