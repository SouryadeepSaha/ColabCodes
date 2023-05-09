//
// Created by Saha on 08/05/2023.
//

#include "GMG_solver.h"

GMG_solver::GMG_solver(int nlvl, int ncycle, double h, double domain, int npre, int npost, int ncoarse) {
    this->poissonGmg = new Poisson_GMG(nlvl, ncycle, h, domain);
    this->npre = npre;
    this->npost = npost;
    this->ncoarse = ncoarse;

}

void GMG_solver::VCycle(int lvl) {

    std::vector<Grid *> sol = this->poissonGmg->get_SOL();
    std::vector<Grid *> f = this->poissonGmg->get_RHS();
    if(lvl == 1){//Coarsest Grid
//        std::cout<<"Lvl Vcycle: "<<lvl<<std::endl;//TODO: check!
        for (int i=0; i<ncoarse;i++){
            this->poissonGmg->GS_Smoothing(lvl);
        }
    }
    else{
//        std::cout<<"Lvl Vcycle: "<<lvl<<std::endl;//TODO: check!
        //pre-smoothing
        for (int i=0; i<npre;i++){
            this->poissonGmg->GS_Smoothing(lvl);
        }

        // compute and restrict the residual
        this->poissonGmg->Restrict_residual(lvl);

        // initialize the coarse solution to zero
        for(unsigned long i=0;i<sol[lvl-1]->get_nrow();i++){
            for(unsigned long j=0;j<sol[lvl-1]->get_ncol();j++){
                sol[lvl-1]->set_val(i,j,0.0);
            }
        }
        //recursive VCycle
        VCycle(lvl-1);

        // interpolate error and correct fine solution
        this->poissonGmg->Prolongate_correct(lvl-1,sol[lvl],sol[lvl-1]);

        //post-smoothing
        for (int i=0; i<npost;i++){
            this->poissonGmg->GS_Smoothing(lvl);
        }

    }
}

std::vector<double> GMG_solver::MultiGrid() {
//    std::cout<<"nlv_MultiGrid(): "<<this->poissonGmg->nlv<<std::endl;//TODO: check!
    std::vector<double> residuals(this->poissonGmg->nlv);
    double res_0 = get_res0();
    residuals[0] = res_0;
    double _res_old{res_0};
    for(int l=1; l<this->poissonGmg->nlv;++l){
        for(int i = 0; i<this->poissonGmg->ncycle; i++){
            VCycle(l);
            _res_old = this->poissonGmg->residual(l);

        }//TODO:Print res
        residuals[l] = _res_old;
        std::cout<<"MultiGrid()_residual: "<<_res_old<<" level: "<< l <<std::endl;
        if(l<this->poissonGmg->nlv-1){
            this->poissonGmg->Prolongation(l-1,this->poissonGmg->get_SOL()[l],this->poissonGmg->get_RHS()[l-1]);
        }
    }


//    std::destroy(this->poissonGmg->get_SOL().begin(), this->poissonGmg->get_SOL().end());
//    std::destroy(this->poissonGmg->get_RHS().begin(), this->poissonGmg->get_RHS().end());
    return residuals;
}

double GMG_solver::get_res0() {
    std::vector<Grid *> f = this->poissonGmg->get_RHS();

    for(int i=this->poissonGmg->nlv-1;i>0 ;i--){
        auto ufine = f[i];
        auto ucoarse = f[i-1];
        this->poissonGmg->Restriction(ufine,ucoarse);
    }
    double res_0 = this->poissonGmg->residual(0);
//    std::cout<<"res_0: "<<res_0<<std::endl;//TODO: check!
    return res_0;
}

void GMG_solver::saveSol(const std::string& filename) {
    /***
     * save Coarsest and Finest grid solution and RHS with appropriate filenames.
     */
    std::ofstream data;
    this->poissonGmg->get_SOL()[0]->saveGrid("./results/Coarsest_"+filename,data);
    this->poissonGmg->get_RHS()[0]->saveGrid("./results/CoarsestF_"+filename,data);
    this->poissonGmg->get_SOL()[this->poissonGmg->nlv-1]->saveGrid("./results/Finest_"+filename,data);
    this->poissonGmg->get_RHS()[this->poissonGmg->nlv-1]->saveGrid("./results/FinestF_"+filename,data);
}