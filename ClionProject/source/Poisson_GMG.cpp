//
// Created by Saha on 07/05/2023.
//

#include "Poisson_GMG.h"


Poisson_GMG::Poisson_GMG(int levels, int iters_V, double h, double domain) {
    nlv = levels;
    ncycle = iters_V;
    invh2 = std::vector<double>(nlv);
    range = domain;
    this->h = h;
//    A = std::vector<Grid*>(nlv,new Grid(1,1));
    sol = std::vector<Grid*>(nlv);//, new Grid(1, 1));
    F = std::vector<Grid*>(nlv);//, new Grid(1, 1));
    for(int i = 0; i<nlv; i++){
        int rlvl = static_cast<int>(pow(2,i)+1);//each level has 2^lvl+1 rows and cols
        int clvl = static_cast<int>(pow(2,i)+1);
        double gridsize = this->h/pow(2,i);
        sol[i] = new Grid(rlvl, clvl);//->resize(rlvl, clvl);
        F[i] = new Grid(rlvl, clvl);//->resize(rlvl, clvl);
        F[i] = ProblemDesc::source_function(F[i], gridsize);//define RHS
        invh2[i] = static_cast<double>(1.0/gridsize/gridsize);//TODO: Check if the 4 needs to be included!
//        std::cout<<"R: "<<rlvl<<" C: "<<clvl<<" Gridsize: "<<gridsize<<std::endl;//TODO: done!
//        F[i]->printGrid();
    }
}
void Poisson_GMG::treatBoundaryNM(int lvl) {
    this->sol[lvl] = ProblemDesc::boundary_conditionsNM(this->sol[lvl], this->invh2[lvl]);

}
void Poisson_GMG::treatBoundaryDR(int lvl){
    this->sol[lvl] = ProblemDesc::boundary_conditionsDR(this->sol[lvl], this->F[lvl], this->invh2[lvl]);
}
void Poisson_GMG::GS_Smoothing(int lvl) {
    double denom = 1.0/(4*this->invh2[lvl]);//diagonal terms of the A matrix

    treatBoundaryDR(lvl);
    //TODO: treatBoundaryNM(lvl);
    for(unsigned long i=1; i< this->sol[lvl]->get_nrow() - 1; i++) {
        for (unsigned long j = 1; j < this->sol[lvl]->get_ncol() - 1; j++) {
            auto val = (this->F[lvl]->get_val(i, j) + \
                               this->invh2[lvl] * (\
                               +this->sol[lvl]->get_val(i - 1, j)\
                               + this->sol[lvl]->get_val(i + 1, j)\
                               + this->sol[lvl]->get_val(i, j - 1)\
                               + this->sol[lvl]->get_val(i, j + 1))\
                               ) * denom;
            this->sol[lvl]->set_val(i, j, val);
        }
    }
    treatBoundaryDR(lvl);
    //TODO: treatBoundaryNM(lvl);
}



#define RES(i,j) (this->F[lvl]->get_val(i,j)/this->invh2[lvl]\
                  +this->sol[lvl]->get_val(i-1,j)\
                  +this->sol[lvl]->get_val(i+1,j)\
                  +this->sol[lvl]->get_val(i,j-1)\
                  +this->sol[lvl]->get_val(i,j+1)\
                  -4.0*this->sol[lvl]->get_val(i,j)) //TODO: f[lvl](i,j)*h^2
RealType Poisson_GMG::residual(int lvl) {
    /***
     * Residual of current level.
     */
    double res(0.0), rf ;
    for(unsigned long i=1; i< this->sol[lvl]->get_nrow() - 1; i++){
        for (unsigned long  j = 1; j < this->sol[lvl]->get_ncol() - 1; j++) {
            rf = this->invh2[lvl]*RES(i,j);//r_lvl = F_lvl-A_lvl*(u_lvl)^iter_k
            res += rf*rf;

        }
    }
    return sqrt(res)/(this->sol[lvl]->get_nrow() * this->sol[lvl]->get_ncol()); //L2 norm
}

void Poisson_GMG::Restriction(Grid *fine, Grid *coarse) {
    /***
     * Interpolate from fine grid to coarse grid.
     */
    // loop over coarse grid points
    for (unsigned long  i = 1; i < coarse->get_nrow()-1; i++) {
        unsigned long  fi = 2*i;
        for( unsigned long j = 1; j< coarse->get_ncol()-1; j++){
            unsigned long fj = 2*j;
            auto IhH = fine->get_val(fi,fj)+fine->get_val(fi-1,fj)\
                       +fine->get_val(fi,fj-1)+fine->get_val(fi-1,fj-1);
            coarse->set_val(i,j,0.25*IhH);
        }

    }
}

void Poisson_GMG::Restrict_residual(int lvl) {
    /***
     * Restrict residual during VCycle from fine to next coarse grid.
     * If lvl_fine = L = lvl (in arguments); lvl_coarse = L-1 .
     */
    // loop over coarse grid points
    for(unsigned long i=1; i< this->F[lvl-1]->get_nrow() - 1; i++){
        unsigned long fi = 2*i;
        for(unsigned long j=1; j< this->F[lvl-1]->get_ncol() - 1; j++){
            unsigned long fj = 2*j;
            auto IhH = this->invh2[lvl]*(RES(fi,fj)+RES(fi-1,fj)+RES(fi,fj-1)+RES(fi-1,fj-1));
            this->F[lvl-1]->set_val(i,j,0.25*IhH);
//            cout<<"IhH_RestrictResidual(): "<<IhH<<endl;
        }

    }
}

void Poisson_GMG::Prolongation(int lvl, Grid *ufine, Grid *ucoarse) {
    /***
     * Interpolate from coarse to fine grid.
     * If lvl_fine = L; lvl_coarse = L-1 = lvl (in arguments).
     */
    // loop over coarse grid points
    for (unsigned long  i = 1; i < ucoarse->get_nrow() - 1; i++) {
        unsigned long fi = 2 * i;
        for (unsigned long j = 1; j < ucoarse->get_ncol() - 1; j++) {
            unsigned long fj = 2 * j;
            auto val = ucoarse->get_val(i, j);
            ufine->set_val(fi, fj, val); //~ufine ( fi , fj ) = v ;
            ufine->set_val(fi - 1, fj, val);
            ufine->set_val(fi, fj - 1, val);
            ufine->set_val(fi - 1, fj - 1, val);
        }
    }
}

void Poisson_GMG::Prolongate_correct(int lvl, Grid *ufine, Grid *ucoarse) {
    /***
     * Interpolate the approximated error from the coarse grid and the correction of the
     * current solution to the fine grid.
     * If lvl_fine = L; lvl_coarse = L-1 = lvl (in arguments).
     */
    for (unsigned long  i = 1; i < ucoarse->get_nrow() - 1; i++) {
        unsigned long fi = 2 * i;
        for (unsigned long j = 1; j < ucoarse->get_ncol() - 1; j++) {
            unsigned long fj = 2 * j;
            auto val = ucoarse->get_val(i, j);
            ufine->set_val(fi, fj, ufine->get_val(fi,fj)+val);//~ufine ( fi , fj ) += val ;
            ufine->set_val(fi - 1, fj, ufine->get_val(fi-1,fj)+val);
            ufine->set_val(fi, fj - 1, ufine->get_val(fi,fj-1)+val);
            ufine->set_val(fi - 1, fj - 1, ufine->get_val(fi-1,fj-1)+val);
        }
    }
}

