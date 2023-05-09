#include "Grid.h"

Grid::Grid(unsigned long rows, unsigned long cols): row(rows),col(cols){
    grid = std::vector<RealType>(row*col,0);
}
RealType Grid::get_val(unsigned long  r, unsigned long  c){
    return grid.at(r+this->col*c);
}
void Grid::set_val(unsigned long  r, unsigned long  c, RealType value){
    grid.at(r+this->col*c) = value;
}

RealType Grid::get_val_T(unsigned long r, unsigned long c){
    return grid.at(r*this->row+c);
}
void Grid::set_val_T(unsigned long r, unsigned long c, RealType value){
    grid.at(r*this->row+c) = value;
}
void Grid::saveGrid(const std::string& filename, std::ofstream& data) {
    data.open(filename, std::ios::out);
    data<<"#X\tY\tZ\n";
    for (unsigned long x = 0; x <this->get_ncol(); x++) {
        for (unsigned long y = 0; y <this->get_nrow(); y++) {
            data << x << "\t"<< y << "\t" << this->get_val(y,x) << "\n";
        }
        data <<"\n";
    }
    data.close();
}

void Grid::printGrid() {
    for (unsigned long x = 0; x <this->get_ncol(); x++) {
        for (unsigned long y = 0; y <this->get_nrow(); y++) {
            cout << x << "\t"<< y << "\t" << this->get_val(y,x) << "\n";
        }
        cout <<"\n";
    }
}