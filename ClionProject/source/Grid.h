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
    RealType get_val(unsigned long r, unsigned long c);
    void set_val(unsigned long  r, unsigned long  c, RealType value);
    RealType get_val_T(unsigned long r, unsigned long c);
    void set_val_T(unsigned long r, unsigned long c, RealType value);
    int get_nrow(){
        return static_cast<unsigned long>(this->grid.size()/col);
    }
    int get_ncol(){
        return static_cast<unsigned long>(this->grid.size()/row);
    }
    void set_row(int r){
        row = r;
        this->grid.resize(row*col);
    }
    void set_col(int c){
        col = c;
        this->grid.resize(row*col);
    }
    std::vector<RealType>* get_grid(){
        return &grid;
    }
//    void resize(int r, int c){
//        if(this->grid.empty()){
//            std::cout<<"Grid not initialised!"<<std::endl;
//            return;
//        }
//        row = r;
//        col = c;
//        this->grid.resize(row*col);
//    }
    void saveGrid(const std::string& filename, std::ofstream& data);
    void printGrid();
private:
    std::vector<RealType> grid;
    unsigned long row;
    unsigned long col;

};


#endif //GRID_H
