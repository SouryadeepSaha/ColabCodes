
#include "Poisson_FD.h"
#include <math.h>
#include <numbers>


Poisson_FD::Poisson_FD(int N, int iterations, RealType range){
    grid_pts = N+1;
    iter = iterations;
    h = (range/static_cast<double>(N));
    f = new Grid(grid_pts,grid_pts);
    u_h = new Grid(grid_pts,grid_pts);
}
void Poisson_FD::source_function() {
    for (int x = 0; x < this->grid_pts; x++) {
        this->f->set_val(x, grid_pts-1, sin(std::numbers::pi * x * this->h));
    }
}

void Poisson_FD::source_function2() {
    for (int y = 0; y < this->grid_pts; y++) {
        for (int x = 0; x < this->grid_pts; x++) {
            this->f->set_val(x, y, 2 * std::numbers::pi * std::numbers::pi *\
                    cos(std::numbers::pi * x * this->h) * cos(std::numbers::pi * y * this->h));

        }
    }
}

void Poisson_FD::boundary_conditions() {
//    g = new Grid(this->grid_pts,this->grid_pts);

    for (int y = 0; y < this->grid_pts; y++) {
        for (int x = 0; x < this->grid_pts; x++) {
            if (y == 0 || y == this->grid_pts - 1) {
                if (x == 0 || x == this->grid_pts - 1) {
                    this->u_h->set_val(x, y, cos(std::numbers::pi * x * this->h) * \
                                        cos(std::numbers::pi * y * this->h));
                }
            }
        }
    }
}

