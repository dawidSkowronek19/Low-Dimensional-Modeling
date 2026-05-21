#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <iomanip>
#include <sstream>

struct Init{
    double L;
    int N;
    //signle
    double a_s;
    double V0_s;
    //gauss
    double sigma;
    double V0_g;
    //double
    double a_d;
    double b_d;
    double V0_d;

    double dx=L/(N-1);
    double t = 1.0/(2.0*dx*dx);
};


Init init0={
    5.0, //L
    500, //N
    1.0, //a_s
    10.0, //V0_s

    0.5, //sigma
    10.0, //V0_g

    1.0, //a_d
    0.5, //b
    10.0, //V0_d


};
const std::complex<double> i = {0.0, 1.0};

void V_init(std::vector<double> &V, std::string potential_type)
{
    if (potential_type=="single")
    {
        for (int x_idx=0; x_idx<init0.N; x_idx++)
        {
            double x=x_idx*init0.dx;
            if (std::abs(init0.L/2-x)<init0.a_s/2)
            {
                V[x_idx]=init0.V0_s;
            }
            else
            {
                V[x_idx]=0.0;
            }
        }
    }

    else if (potential_type=="gaussian")
    {
        for (int x_idx=0; x_idx<init0.N; x_idx++)
        {
            double x = x_idx*init0.dx;
            V[x_idx]=init0.V0_g*exp(-(x-init0.L/2)*(x-init0.L/2)/(init0.sigma*init0.sigma));
        }
    }

    else if (potential_type=="double")
    {
        double d = init0.b_d/2+init0.a_d/2;
        for (int x_idx=0; x_idx<init0.N; x_idx++)
        {
            double x=x_idx*init0.dx;
            if (std::abs(init0.L/2-d-x)<init0.a_d/2.0 || std::abs(init0.L/2+d-x)<init0.a_d/2.0)
            {
                V[x_idx]=init0.V0_d;
            }
            else
            {
                V[x_idx]=0.0;
            }
        }
    }

}

std::vector<std::complex<double>> thomas(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b, std::vector<std::complex<double>> c, std::vector<std::complex<double>> d)
{
    /*
        Ax=d => a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i  
    */
    
    int n = init0.N;
    std::vector<std::complex<double>> x(n);

    // Eliminacja w przód
    c[0] /= b[0];
    d[0] /= b[0];

    for (int I = 1; I < n; I++) {
        std::complex<double> m = b[I] - a[I] * c[I-1];
        if (I < n - 1) c[I] /= m;
        d[I] = (d[I] - a[I] * d[I-1]) / m;
    }

    // Podstawienie wsteczne
    x[n-1] = d[n-1];
    for (int I = n - 2; I >= 0; I--) {
        x[I] = d[I] - c[I] * x[I+1];
    }

    return x;
}



int main()
{

    std::vector<double> V(init0.N);
    V_init(V, "single");

    std::vector<std::complex<double>> diag(init0.N,0.0);
    std::vector<std::complex<double>> e(init0.N, 0.0);
    std::vector<std::complex<double>> S(init0.N, 0.0);
    std::vector<std::complex<double>> psi(init0.N, 0.0);

    double E=0.05;
    double E_max=50.0;
    double dE=0.05;


    std::ofstream fileE("E.dat");
    while (E<E_max)
    {

        double k=sqrt(2.0*E);
        std::complex<double> s=0;
        std::complex<double> SIGMA=0.0;

        for (int x_idx=0; x_idx<init0.N; x_idx++)
        {
            if (x_idx==0)
            {
                s=init0.t*(1.0-exp(2.0*i*k*init0.dx));
                SIGMA=-init0.t*exp(i*k*init0.dx);
            }
            else if(x_idx==init0.N-1)
            {
                s=0.0;
                SIGMA=-init0.t*exp(i*k*init0.dx);
            }
            else
            {
                s=0.0;
                SIGMA=0.0;
            }

            S[x_idx]=s;
            diag[x_idx]=SIGMA+2.0*init0.t + V[x_idx] - E;
            e[x_idx]=-init0.t;
        }

        psi=thomas(e,diag,e,S);

        double T=std::abs(psi[init0.N-1])*std::abs(psi[init0.N-1]);
        fileE<<E<<" "<<T<<" "<<"\n";

        if (T>0.98)
        {
            std::ostringstream filename;
            filename << "./outdir/PSI_E" << std::fixed << std::setprecision(4)<< E << ".dat";

            std::ofstream file(filename.str());

            for (int x_idx=0; x_idx<init0.N; x_idx++)
            {
                file<<x_idx*init0.dx<<" "<<std::abs(psi[x_idx])<<"\n";
            }

            file.close();
        }

        E+=dE;
    }
    fileE.close();
    
    return 0;
}