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
#include <fstream>
#include <iostream>

using namespace std;
class Grid {

public:
    Grid(unsigned long rows, unsigned long cols);
    Grid(unsigned long rows, unsigned long cols, RealType value);
    RealType get_val(unsigned long r, unsigned long c);
    void set_val(unsigned long  r, unsigned long  c, RealType value);
    unsigned long get_nrow(){
        return row;//static_cast<unsigned long>(this->grid.size()/col);
    }
    unsigned long get_ncol(){
        return col;//static_cast<unsigned long>(this->grid.size()/row);
    }
    double L2norm();

//    void set_row(int r){
//        row = r;
//        grid.resize(row*col);
//    }
//    void set_col(int c){
//        col = c;
//        grid.resize(row*col);
//    }

    std::vector<RealType>* get_grid(){
        return &grid;
    }
//    void resize(int r, int c){
//        if(this->grid.empty()){
//            std::cout<<"Grid not initialised!"<<std::endl;
//            return;
//        }
//        this->row = r;
//        this->col = c;
//        this->grid.resize(row*col);
//    }
    void saveGrid(const std::string& filename, std::ofstream& data);
    void printGrid();
private:
    std::vector<RealType> grid;
    // RealType *grid;
    unsigned long row;
    unsigned long col;

};

inline RealType Grid::get_val(unsigned long  r, unsigned long  c){
    //assert(r<row);
    //return grid[c+col*r];
    return grid.at(r*this->col+c);
}

inline void Grid::set_val(unsigned long  r, unsigned long  c, RealType value){
    grid.at(r*this->col+c) = value;
}


#endif //GRID_H
