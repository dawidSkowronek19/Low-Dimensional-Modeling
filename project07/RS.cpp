#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <complex>
#include "Parameters.h"
//======================= PARAMETERS AND STRUCTURES ==============
const int M=4;
struct SparseMatrix{
    Eigen::SparseMatrix<std::complex<double>> Matrix;
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
};

const std::complex<double> I={0,1};

const double dx=init.a/(init.N_x);
const double dk=(2.0*M_PI/init.a)/(init.N_k-1);
const double hbar=1.0;
const double m=1.0;
//================== POTENTIAL ==================
void V_init(std::vector<double> &V)
{
    std::cout<<"COMPUTING POTENTIAL START\n";
    for (int x_idx=0; x_idx<init.N_x; x_idx++)
    {
        double x = -init.a/2+x_idx*dx;
        
        if (std::abs(x)<=init.w/2)
        {
            V[x_idx]=-init.V0;
        }
        else
        {
            V[x_idx]=0.0;
        }
    }
    std::cout<<"COMPUTING POTENTIAL DONE\n";
}
//======================== SOLVER FOR GIVEN K ===============
Eigen::VectorXd E_kStep(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> &solver, const std::vector<double> &V, SparseMatrix &H, double k)
{
    double phi=k*dx;
    double t=hbar*hbar/(2.0*m*dx*dx);

    double D_cst = t*(2.0+phi*phi);
    std::complex<double> L_cst= t*(I*phi-1.0);
    std::complex<double> U_cst=t*(-I*phi-1.0);
    std::cout<<"Building H_matrix start\n";
    for(int i=0; i<init.N_x; ++i) 
    {
        H.triplets.push_back({i,i, D_cst+V[i]});
        if(i < init.N_x-1) 
        {
            H.triplets.push_back({i,i+1, U_cst});
            H.triplets.push_back({i+1, i, L_cst});
        }
    }

    H.triplets.push_back({0, init.N_x-1, L_cst});
    H.triplets.push_back({init.N_x-1,0, U_cst});
    H.Matrix.setFromTriplets(H.triplets.begin(), H.triplets.end());

    std::cout<<"Building H_matrix done\n";

    std::cout<<"Diagonalizing matrix H(k)\n";
    solver.compute(H.Matrix, Eigen::EigenvaluesOnly);
    std::cout<<"Diagonalization done\n";
    return solver.eigenvalues();

}

int main()
{

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    std::vector<double> E(init.N_x*init.N_k,0.0), V(init.N_x, 0.0);
    V_init(V);
    Eigen::VectorXd epsilon_k;
    SparseMatrix H;
    H.Matrix.resize(init.N_x, init.N_x);
    

    std::cout<<"SOLVER START\n\n";
    for (int k_idx=0; k_idx<init.N_k; k_idx++)
    {
        std::cout<<"STEP k_idx="<<k_idx<<"\n";
        double k=-M_PI/init.a + k_idx*dk;
        epsilon_k=E_kStep(solver, V, H, k);

        std::cout<<"Saving eigenvalues start\n";
        for (int x_idx=0; x_idx<init.N_x; x_idx++)
        {
            int l = x_idx + init.N_x*k_idx;
            E[l]=epsilon_k(x_idx);
        }
        std::cout<<"Saving eigenvalues done\n";

        H.triplets.clear();

    }
    std::cout<<"SOLVER DONE\n\n";


    std::cout<<"Saving to file start\n";
    std::ofstream file("bands.dat");
    for (int k_idx=0; k_idx<init.N_k; k_idx++)
    {
        double k=-M_PI/init.a +k_idx*dk;
        file<<k<<" ";
        for (int x_idx=0; x_idx<M; x_idx++)
        {
            int l = x_idx+init.N_x*k_idx;
            file<<E[l]<<" ";
        }
        file<<"\n";
    }
    std::cout<<"Saving to file done\n";
    return 0;
}