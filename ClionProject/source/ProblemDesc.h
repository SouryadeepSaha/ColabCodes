//
// Created by Saha on 08/05/2023.
//

#ifndef PROBLEMDESC_H
#define PROBLEMDESC_H
#include "Grid.h"
#include <cmath>
#include <numbers>

class ProblemDesc {
public:
    ProblemDesc();
    static Grid* source_function(Grid *f, double h);
    static Grid* source_function2(Grid *f, double h);
    static Grid* boundary_conditions(Grid *u_h, double h) ;
    static Grid* boundary_conditionsDR(Grid *u_h,Grid *g);
    static Grid* boundary_conditionsNM(Grid *u_h, double h);
    static Grid* exactSolution(unsigned long nrows, unsigned long ncols, double h);


};


#endif //PROBLEMDESC_H
