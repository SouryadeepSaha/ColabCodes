
#include "dependencies.h"
#include "chrono"
using namespace std;
int main() {
    cout << "Hello, World!" << endl;
//    int N{100}, iter{10000};
//    double range {1.0};
//    auto poissonGS = new Poisson_FD(N, iter, range);
//    poissonGS->source_function();
//    poissonGS->boundary_conditions();
    auto start = chrono::high_resolution_clock::now();
//    poissonGS->gauss_seidel();
    auto stop = chrono::high_resolution_clock::now();
//    auto duration = duration_cast<chrono::seconds>(stop - start);
//    cout<<"Gauss Seidel execution time for "<<N<<" grid points and "<<iter<<\
//        " iterations is "<<duration.count()<<" seconds."<<endl;
//
//    string filename = "./results/solution.txt";
//    ofstream data;
//
//    data.open(filename, ios::out);
//    data<<"#X\tY\tZ\n";
//    for (int x = 0; x <=N; x++) {
//        for (int y = 0; y <=N; y++) {
//            data << x << "\t"<< y << "\t" << poissonGS->get_Uh()->get_val(y,x) << "\n";
//        }
//        data <<endl;
//    }
//    data.close();
//
//    auto poissonGS2 = new Poisson_FD(N, iter, range);
//    poissonGS2->source_function2();
//    poissonGS2->boundary_conditions2();
//    poissonGS2->gauss_seidel();
//
//    filename = "./results/solution2.txt";
////    ofstream data;
//
//    data.open(filename, ios::out);
//    data<<"#X\tY\tZ\n";
//    for (int x = 0; x <=N; x++) {
//        for (int y = 0; y <=N; y++) {
//            data << x << "\t"<< y << "\t" << poissonGS2->get_Uh()->get_val(y,x) << "\n";
//        }
//        data <<endl;
//    }
//    data.close();

    int Nlvl{8}, Vcycle{5};
    start = chrono::high_resolution_clock::now();

    auto mgsolve = new GMG_solver(Nlvl,Vcycle,1.0,1.0,5,10,10);
    auto residuals = mgsolve->MultiGrid();

    stop = chrono::high_resolution_clock::now();
    auto duration2 = duration_cast<chrono::seconds>(stop - start);
    cout<<"Multigrid execution time for "<<Nlvl<<" levels and "<<Vcycle<<\
        " Vcycles is "<<duration2.count()<<" seconds."<<endl;

    mgsolve->saveSol("MGsolutionsTest.txt");

    for(auto i=0; i<residuals.size();i++){
       cout<<"Residual for level "<<i<<" : "<<residuals.at(i)<<endl;
    }
//    cout<<"Starting residual: "<<residuals.at(1)<<endl;
//    cout<<"Final residual: "<<residuals.at(residuals.size()-1)<<endl;

    return 0;
}

