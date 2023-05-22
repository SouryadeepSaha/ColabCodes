//
// Created by Saha on 07/05/2023.
//

#include "Poisson_GMG.h"


Poisson_GMG::Poisson_GMG(int levels, int iters_V, double h, double domain) {
    nlv = levels;
    ncycle = iters_V;
    h2 = std::vector<double>(nlv);
    range = domain;
    this->h = h; //h;
    U_exact = std::vector<Grid*>(nlv,new Grid(1,1));
    sol = std::vector<Grid*>(nlv);//, new Grid(1, 1));
    F = std::vector<Grid*>(nlv);
    Res = std::vector<Grid*>(nlv);//, new Grid(1, 1));
    for(int i = 0; i<nlv; i++){
        int rlvl = static_cast<int>(pow(2,i)+1);//each level has 2^lvl+1 rows and cols
        int clvl = static_cast<int>(pow(2,i)+1);
        double gridsize = this->h/pow(2,i);
        //setup all the grids
        sol[i] = new Grid(rlvl, clvl,1.0);//->resize(rlvl, clvl);
        F[i] = new Grid(rlvl, clvl);//->resize(rlvl, clvl);
        U_exact[i] =  new Grid(rlvl, clvl);
        Res[i] = new Grid(rlvl, clvl);
        //define RHS and exact solution
//        F[i] = ProblemDesc::source_function2(F[i], gridsize);
//        U_exact[i] = ProblemDesc::exactSolution(rlvl,clvl,gridsize);
        h2[i] = static_cast<double>(gridsize * gridsize);//

    }
}
void Poisson_GMG::treatBoundaryNM(int lvl) {
    this->sol[lvl] = ProblemDesc::boundary_conditionsNM(this->sol[lvl], this->h2[lvl]);

}
void Poisson_GMG::treatBoundaryDR(int lvl){
    this->sol[lvl] = ProblemDesc::boundary_conditionsDR(this->sol[lvl], this->F[lvl]);
}
void Poisson_GMG::GS_Smoothing(int lvl) {

    treatBoundaryDR(lvl);
    //TODO: treatBoundaryNM(lvl);
    for(unsigned long i=1; i< this->sol[lvl]->get_nrow() - 1; i++) {
        for (unsigned long j = 1; j < this->sol[lvl]->get_ncol() - 1; j++) {
            auto val = (this->h2[lvl]*this->F[lvl]->get_val(i, j) + \
                               (this->sol[lvl]->get_val(i - 1, j)\
                               + this->sol[lvl]->get_val(i + 1, j)\
                               + this->sol[lvl]->get_val(i, j - 1)\
                               + this->sol[lvl]->get_val(i, j + 1))\
                               ) /4.0;
            this->sol[lvl]->set_val(i, j, val);
        }
    }
    treatBoundaryDR(lvl);
    //TODO: treatBoundaryNM(lvl);
}



#define RES(i,j) (F[lvl]->get_val(i,j)\
                  +(sol[lvl]->get_val(i-1,j)\
                  +sol[lvl]->get_val(i+1,j)\
                  +sol[lvl]->get_val(i,j-1)\
                  +sol[lvl]->get_val(i,j+1)\
                  -4.0*sol[lvl]->get_val(i,j))/h2[lvl])//F_lvl-A_lvl*(u_lvl)

void Poisson_GMG::residual(int lvl) {
    /***
     * Storing residual of current level.
     */
    //double res(0.0), rf(0.0) ;
    for(unsigned long i=1; i< sol[lvl]->get_nrow() - 1; i++){
        for (unsigned long  j = 1; j < sol[lvl]->get_ncol() - 1; j++) {

            Res[lvl]->set_val(i,j,
                              F[lvl]->get_val(i,j)
                               + (sol[lvl]->get_val(i-1,j)
                               + sol[lvl]->get_val(i+1,j)
                               + sol[lvl]->get_val(i,j-1)
                               + sol[lvl]->get_val(i,j+1)
                               -4.0*sol[lvl]->get_val(i,j)) / h2[lvl]);
            //r_lvl = F_lvl-A_lvl*(u_lvl)^iter_k
        }
    }
    //return sqrt(res)/(this->sol[lvl]->get_nrow() * this->sol[lvl]->get_ncol()); //L2 norm
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
            auto IhH = fine->get_val(fi,fj+1)+fine->get_val(fi+1,fj)\
                       +fine->get_val(fi,fj-1)+fine->get_val(fi-1,fj);
            coarse->set_val(i,j,0.25*IhH);
        }

    }
}

void Poisson_GMG::Restrict_residual(int lvl) {
    /***
     * Restrict residual during VCycle from fine to next coarse grid.
     * lvl (in arguments) = lvl_fine = L; lvl_coarse = L-1 .
     */
    Restriction(Res[lvl],Res[lvl-1]);
    // loop over coarse grid points
    for(unsigned long i=1; i< this->F[lvl-1]->get_nrow() - 1; i++){
        unsigned long fi = 2*i;
        for(unsigned long j=1; j< this->F[lvl-1]->get_ncol() - 1; j++){
            unsigned long fj = 2*j;
            auto IhH = (RES(fi, fj+1) + RES(fi + 1, fj) + RES(fi, fj - 1) + RES(fi - 1, fj ));
            this->F[lvl-1]->set_val(i,j,0.25*IhH);
//            cout<<"IhH_RestrictResidual(): "<<IhH<<endl;
        }

    }
}

void Poisson_GMG::Prolongation(int lvl, Grid *ufine, Grid *ucoarse) {
    /***
     * Interpolate from coarse to fine grid.
     * If lvl_fine = L; lvl (in arguments) = lvl_coarse = L-1.
     */
    // loop over coarse grid points
    for (unsigned long  i = 1; i < ucoarse->get_nrow() - 1; i++) {
        unsigned long fi = 2 * i;
        for (unsigned long j = 1; j < ucoarse->get_ncol() - 1; j++) {
            unsigned long fj = 2 * j;
            auto val = ucoarse->get_val(i, j);
            ufine->set_val(fi, fj+1, val); //~ufine ( fi , fj ) = v ;
            ufine->set_val(fi + 1, fj, val);
            ufine->set_val(fi, fj - 1, val);
            ufine->set_val(fi - 1, fj, val);
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
            cout<<"ProlongateCorrect()_Coarse Grid value: "<<val<<" level: "<<lvl<<endl;
            ufine->set_val(fi, fj+1, ufine->get_val(fi,fj)+val);//~ufine ( fi , fj ) += val ;
            ufine->set_val(fi + 1, fj, ufine->get_val(fi-1,fj)+val);
            ufine->set_val(fi, fj - 1, ufine->get_val(fi,fj-1)+val);
            ufine->set_val(fi - 1, fj, ufine->get_val(fi-1,fj-1)+val);
        }
    }
}

double Poisson_GMG::L2norm(Grid *grid1, Grid *grid2) {
    if(grid1->get_nrow() == grid2->get_nrow() && grid1->get_ncol() == grid2->get_ncol()){
        unsigned long row = grid1->get_nrow();
        unsigned long col = grid1->get_ncol();
        double norm2{0.0};
        for (unsigned long i = 0; i < row; ++i) {
            for (unsigned long j = 0; j < col; ++j) {
                 double val = grid1->get_val(i,j)-grid2->get_val(i,j);
                 norm2 +=val*val;
            }
        }
        return sqrt(norm2);
    }
    else
        return -1.0;

}
