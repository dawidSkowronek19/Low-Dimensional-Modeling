#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <complex>
#include "./Parameters.h"
const Parameters init={
    1,//a
    100,//V0
    200, //N_k
    1001, //N_x
};

const int p=1;
double x_center=0.0;

int M=3; //nb of states
const std::string dis_usage="High";

double x_max=x_center+(p+1.5)*init.a;
double x_min=x_center-(p+1.5)*init.a;

const double dx=(x_max-x_min)/(init.N_x-1);
const double dk=(2.0*M_PI/init.a)/(init.N_k-1);
const double hbar=1.0;
const double m=1.0;
const std::complex<double> I ={0.0, 1.0};

double V_local(double x)
{
    if (std::abs(x-x_center) <= init.w/2.0) 
    {
        return -init.V0;
    }
    return 0.0;
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> solveLocalSchrodinger()
{
    std::cout<<"Local Schrodinger solver...\n";
    Eigen::VectorXd diag(init.N_x);
    Eigen::VectorXd subdiag(init.N_x-1);

    subdiag.setConstant(-hbar*hbar/(2.0*m*dx*dx));

    for (int i=0; i<init.N_x; i++)
    {
        double x=x_min+i*dx;
        diag(i)=2.0*hbar*hbar/(2.0*m*dx*dx) + V_local(x);
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    solver.computeFromTridiagonal(diag, subdiag, Eigen::ComputeEigenvectors);

    if (solver.info() != Eigen::Success) 
    {
        std::cout << solver.info()<<"\n";
    }

    std::cout<<"Local Schrodinger solver done\n";
    return {solver.eigenvectors().leftCols(M), solver.eigenvalues().head(M)};
}


double V_total(double x)
{
    double V=0.0;
    for (int p_idx=-p; p_idx<=p; p_idx++)
    {
        
        V+=V_local(x-p_idx*init.a);
    }

    return V;
}

std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> HS_matrix(double k, Eigen::MatrixXd &Psi, Eigen::VectorXd &E)
{
    std::complex<double> s=0.0, h=0.0;

    Eigen::MatrixXcd S(M,M), H(M,M);
    S.setConstant(0.0);
    H.setConstant(0.0);



    for (int m=0; m<M; m++)
    {
        for (int n=0; n<M; n++)
        {
            for (int p_idx=-p; p_idx<=p; p_idx++)
            {
                for (int x_idx=0; x_idx<init.N_x; x_idx++)  //x_p = x_min +x_p_idx * dx
                {

                    double x = x_min + x_idx*dx ;
                    double x_p = x - p_idx*init.a;
                    int x_p_idx = std::round((x_p-x_min)/dx);
                    if (x_p_idx >= 0 && x_p_idx < init.N_x) 
                    {
                        s += std::conj(Psi(x_idx, m)) * Psi(x_p_idx, n);
                        h += std::conj(Psi(x_idx, m)) * (E(n) + V_total(x) - V_local(x_p)) * Psi(x_p_idx, n);
                    }
                }

            double dR=p_idx*init.a;
            S(m,n)+=s*exp(-I*k*dR)*dx;
            H(m,n)+=h*exp(-I*k*dR)*dx;

            h=0.0; 
            s=0.0;
            }
        }
    }
    
    return {H, S};

}




int main()
{

    auto[Psi, E] = solveLocalSchrodinger();

    for (int m=0; m<M; m++)
    {
        std::ofstream file_psi("./outdir/Psi_local"+std::to_string(m)+".dat");

        for (int x_idx=0; x_idx<init.N_x; x_idx++)
        {
            double x = x_min+x_idx*dx;
            file_psi<<x<<" "<<Psi(x_idx, m)<<'\n';
        }

        file_psi.close();
    }


    double k_min=-M_PI/init.a;
    double k_max=M_PI/init.a;

    Eigen::MatrixXcd Epsilon(M,init.N_k);
    std::ofstream file_epsilon("./epsilon.dat");
    while(k_min<=k_max)
    {
        std::cout<<"Solver for k= "<<k_min<<"\n";
        auto [H, S] = HS_matrix(k_min, Psi, E);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H, S);
        if(solver.info() != Eigen::Success) 
        {
            std::cout << "Solver failed\n";
            return 1;
        }

        file_epsilon<<k_min<<" ";
        for (int i=0; i<M; i++)
        {
            file_epsilon<<solver.eigenvalues()(i)<<" ";
        }
        file_epsilon<<"\n";
        k_min+=dk;
        H.setConstant(0.0);
        S.setConstant(0.0);
    }
    file_epsilon.close();
    return 0;
}