//
// Created by soury on 01/06/2022.
//

#ifndef GRID_H
#define GRID_H
#ifdef USE_FLOATS
    typedef float RealType;
#else
    typedef double RealType;
#endif
#include <vector>
#include <cassert>

class Grid {

public:
    Grid(int rows, int cols);
    RealType get_val(int r, int c);
    void set_val(int r, int c, RealType value);
    RealType get_val_T(int r, int c);
    void set_val_T(int r, int c, RealType value);
private:
    std::vector<RealType> grid;
    int row;
    int col;

};


#endif //GRID_H
