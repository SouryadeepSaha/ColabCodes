#include "Grid.h"

Grid::Grid(unsigned long rows, unsigned long cols): row(rows),col(cols){
    grid = std::vector<RealType>(row*col,0.0);

    //  grid = new RealType[row*col];
}
Grid::Grid(unsigned long rows, unsigned long cols, RealType value): row(rows),col(cols){
    grid = std::vector<RealType>(row*col,value);
}

double Grid::L2norm() {
    /***
     * Return L2 norm squared for the grid.
     */
    RealType norm2 = 0.0;
    for (unsigned long i = 0; i < row; ++i) {
        for (unsigned long j = 0; j < col; ++j) {
            double val = this->get_val(i,j);
            norm2 +=val*val;
        }
    }
    return norm2;
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