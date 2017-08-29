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

        // // Instead of staring the border coordinates, we need only store the midpoint where the node divides.
        // // This means we have fewer checks that can be done, but we save in time and space.
        // unsigned int mMidX;
        // unsigned int mMidY;

        // Zone de délimitation représentant les limites de ce Tree

        // Vecteur pour stocker les cellules selectionnées dans le Tree
        std::vector<Cell*>* cells;

        // // Constante arbitraire indiquant combien d'éléments peuvent être stockés au maximum dans un noeud
        // static constexpr int QT_NODE_CAPACITY = 15;
    protected:
        Rectangle* rectangle;
};

#endif
