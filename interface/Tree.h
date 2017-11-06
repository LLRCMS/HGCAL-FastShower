#ifndef TREE_H_
#define TREE_H_

#include "Cell.h"
#include "Rectangle.h"

#include <vector>

class Tree {

    public:
        Tree(Rectangle*, int);
        ~Tree();

        bool empty();
        void addCell(Cell*);
        int countCells();
        std::vector<Cell*>* getCells();
        Tree* getLeaf(Point*);
        Tree* getLeaf(float, float);

    private:
        //4 childrens
        Tree* nw;
        Tree* ne;
        Tree* sw;
        Tree* se;

        void subdivide(int);

        // Cells selected in the tree
        std::vector<Cell*>* cells;

    protected:
        Rectangle* rectangle;
};

#endif
