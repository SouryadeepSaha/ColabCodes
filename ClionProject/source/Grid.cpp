#include "Grid.h"

Grid::Grid(int rows, int cols): row(rows),col(cols){

    grid = std::vector<RealType>(row*col,0);
}
RealType Grid::get_val(int r, int c){
    return grid.at(r+this->col*c);
}
void Grid::set_val(int r, int c, RealType value){
    grid.at(r+this->col*c) = value;
}

RealType Grid::get_val_T(int r, int c){
    return grid.at(r*this->col+c);
}
void Grid::set_val_T(int r, int c, RealType value){
    grid.at(r*this->col+c) = value;
}