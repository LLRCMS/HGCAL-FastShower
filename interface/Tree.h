#ifndef __HGCalSimulation_FastShower_Tree_h__
#define __HGCalSimulation_FastShower_Tree_h__

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
    Tree(const Rectangle&, int);
    ~Tree();

    bool empty();
    void addCell(const Cell*);
    int countCells();
    const std::vector<const Cell*> getCells();
    Tree* getLeaf(const Point);
    Tree* getLeaf(float, float);

  private:
    //4 childrens
    std::unique_ptr<Tree> nw_;
    std::unique_ptr<Tree> ne_;
    std::unique_ptr<Tree> sw_;
    std::unique_ptr<Tree> se_;

    void subdivide(int);

    // Cells selected in the tree
    std::vector<const Cell*> cells_;
    Rectangle rectangle_;
};

#endif