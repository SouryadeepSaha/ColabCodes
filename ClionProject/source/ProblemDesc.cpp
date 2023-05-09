//
// Created by Saha on 08/05/2023.
//

#include "ProblemDesc.h"

Grid* ProblemDesc::source_function(Grid *f, double h) {
    unsigned long row = f->get_nrow();
    unsigned long col = f->get_ncol();

    for (unsigned long x = 0; x < col; x++) {
        f->set_val(row-1, x, sin(std::numbers::pi * x * h));
    }
    return f;
}
Grid* ProblemDesc::source_function2(Grid *f, double h) {
    unsigned long row = f->get_nrow();
    unsigned long col = f->get_ncol();
    for (unsigned long y = 0; y < row; y++) {
        for (unsigned long x = 0; x < col; x++) {
            f->set_val(y, x, 2 * std::numbers::pi * std::numbers::pi *\
                    cos(std::numbers::pi * x * h) * cos(std::numbers::pi * y * h));

        }
    }
    return f;
}

Grid* ProblemDesc::boundary_conditionsDR(Grid *u_h, Grid *g,double h) {
    unsigned long row = u_h->get_nrow();
    unsigned long col = u_h->get_ncol();
    for (unsigned long y = 0; y < row; y++) {
        for (unsigned long x = 0; x < col; x++) {
            if (y == 0 || y == row - 1 || x == 0 || x == col - 1) {
                u_h->set_val(y, x, g->get_val(y,x));
//                u_h->set_val(y, x, 0.5*g->get_val(y,x)/numbers::pi/numbers::pi);
            }
        }

    }
    return u_h;
}
Grid* ProblemDesc::boundary_conditionsNM(Grid *u_h, double h) {
    unsigned long row = u_h->get_nrow();
    unsigned long col = u_h->get_ncol();
    //left and right boundary:
    for (unsigned long y = 0; y < row; y++) {
        u_h->set_val(y,0,u_h->get_val(y,1));
        u_h->set_val(y,col-1,u_h->get_val(y,col-2));
    }
    //top and bottom boundary:
    for (unsigned long x = 0; x < col; x++) {
        u_h->set_val(0,x,u_h->get_val(1,x));
        u_h->set_val(row-1,x,u_h->get_val(row-2,x));
    }
    return u_h;
}