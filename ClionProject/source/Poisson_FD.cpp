
#include "Poisson_FD.h"



Poisson_FD::Poisson_FD(unsigned long N, int iterations, RealType range){
    grid_pts = N+1;
    iter = iterations;
    h = (range/static_cast<RealType>(N));
    f = new Grid(grid_pts,grid_pts);
    u_h = new Grid(grid_pts,grid_pts);
}
void Poisson_FD::source_function() {
    for (unsigned long x = 0; x < this->grid_pts; x++) {
        this->f->set_val(grid_pts-1, x, sin(std::numbers::pi * x * this->h));
    }
}

void Poisson_FD::source_function2() {
    for (unsigned long y = 0; y < this->grid_pts; y++) {
        for (unsigned long x = 0; x < this->grid_pts; x++) {
            this->f->set_val(y, x, 2 * std::numbers::pi * std::numbers::pi *\
                    cos(std::numbers::pi * x * this->h) * cos(std::numbers::pi * y * this->h));

        }
    }
}
void Poisson_FD::boundary_conditions() {
    for (int y = 0; y < this->grid_pts; y++) {
        for (int x = 0; x < this->grid_pts; x++) {
            if (y == 0 || y == this->grid_pts - 1 || x == 0 || x == this->grid_pts - 1) {
                this->u_h->set_val(y, x, this->f->get_val(y,x));
            }
        }
    }
}
void Poisson_FD::boundary_conditions2() {
//    g = new Grid(this->grid_pts,this->grid_pts);

    for (unsigned long y = 0; y < this->grid_pts; y++) {
        for (unsigned long x = 0; x < this->grid_pts; x++) {
            if (y == 0 || y == this->grid_pts - 1 || x == 0 || x == this->grid_pts - 1) {
                this->u_h->set_val(y, x, cos(std::numbers::pi * x * this->h) * \
                                        cos(std::numbers::pi * y * this->h));
            }
        }
    }
}

void Poisson_FD::gauss_seidel() {
    for (int i = 0; i < this->iter; ++i) {
        for (unsigned long y = 1; y < this->grid_pts - 1; y++) {
            for (unsigned long x = 1; x < this->grid_pts - 1; x++) {
                RealType update = this->f->get_val(y,x)*this->h*this->h +\
                                  this->u_h->get_val(y-1, x) +\
                                  this->u_h->get_val(y+1, x) +\
                                  this->u_h->get_val(y, x-1) +\
                                  this->u_h->get_val(y, x+1) ;
                this->u_h->set_val(y,x,static_cast<RealType>(update/4.0));
            }
        }
    }
    //return u_h;
}

Grid* Poisson_FD::get_Uh() {
    return this->u_h;
}

Grid* Poisson_FD::get_F() {
    return this->f;
}