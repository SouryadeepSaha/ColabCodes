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
    void residual(int lvl);
    void Restrict_residual(int lvl);
    void Prolongation(int lvl, Grid *ufine, Grid *ucoarse);
    void Prolongate_correct(int lvl, Grid *ufine, Grid *ucoarse);
    double L2norm(Grid* grid1, Grid* grid2);
    double L2residual(int lvl){
        return sqrt(Res[lvl]->L2norm()/(Res[lvl]->get_nrow()*Res[lvl]->get_ncol()));
    }

    std::vector<Grid*> get_SOL(){return sol;}
    std::vector<Grid*> get_RHS(){return F;}
    std::vector<Grid*> get_Exact(){return U_exact;}
    int nlv,ncycle;
    std::vector<double> h2;
    double range;


private:
    std::vector<Grid*> U_exact;
    std::vector<Grid*> sol;
    std::vector<Grid*> F;
    std::vector<Grid*> Res;
    double h;



};


#endif //POISSON_GMG_H
