//
// Created by soury on 31/05/2022.
//

#ifndef POISSON_FD_H
#define POISSON_FD_H
#include "Grid.h"

class Poisson_FD {


public:
    Poisson_FD(int N, int iterations, RealType range);
    void source_function();
    void source_function2();
    void boundary_conditions();
    Grid* gauss_seidel(){return U;}
private:
    int grid_pts{};
    double h{};
    int iter;
    Grid* f;
    Grid* u_h;
    Grid* g;

};


#endif // POISSON_FD_H
