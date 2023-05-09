//
// Created by Saha on 07/05/2023.
//

#ifndef POISSON_GMG_H
#define POISSON_GMG_H
#include "Grid.h"
#include "ProblemDesc.h"

class Poisson_GMG {
public:
    Poisson_GMG(int levels, int iters_V, double h=1.0, double domain=1.0);
    void GS_Smoothing(int lvl);
    void treatBoundaryDR(int lvl);
    void treatBoundaryNM(int lvl);
    void Restriction(Grid *fine, Grid *coarse);
    RealType residual(int lvl);
    void Restrict_residual(int lvl);
    void Prolongation(int lvl, Grid *ufine, Grid *ucoarse);
    void Prolongate_correct(int lvl, Grid *ufine, Grid *ucoarse);
    std::vector<Grid*> get_SOL(){return sol;}
    std::vector<Grid*> get_RHS(){return F;}
    int nlv,ncycle;
    std::vector<double> invh2;
    double range, h;


private:
//    std::vector<Grid*> A;
    std::vector<Grid*> sol;
    std::vector<Grid*> F;



};


#endif //POISSON_GMG_H
