//
// Created by Saha on 08/05/2023.
//

#ifndef GMG_SOLVER_H
#define GMG_SOLVER_H
#include "Poisson_GMG.h"

class GMG_solver {
public:
    GMG_solver(int nlvl, int ncycle, double h, double domain,\
                int npre, int npost, int ncoarse);
    std::vector<double> MultiGrid();
    void saveSol(const std::string& filename);
    int npre{1}, npost{1}, ncoarse{1};
//    bool choice;//TODO!

private:
    double get_res0();
    void VCycle(int lvl);
    Poisson_GMG* poissonGmg;
};


#endif //GMG_SOLVER_H
