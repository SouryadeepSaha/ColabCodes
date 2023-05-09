//
// Created by soury on 31/05/2022.
//

#ifndef POISSON_FD_H
#define POISSON_FD_H
#include "Grid.h"
#include <cmath>
#include <numbers>

class Poisson_FD {


public:
    Poisson_FD(unsigned long N, int iterations, RealType range);
    void source_function();
    void source_function2();
    void boundary_conditions();
    void boundary_conditions2();
    void gauss_seidel();
    Grid* get_Uh();
    Grid* get_F();
private:
    unsigned long grid_pts{};
    double h{};
    int iter;
    Grid* f;
    Grid* u_h;
//    Grid* g;

};


#endif // POISSON_FD_H
