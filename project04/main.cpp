#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>

double h=1.0, m=1.0;
const std::complex<double> i ={0.0, 1.0};

struct Init{

    double dx;
    double dt;
    double L;
    double V0;
    double k0;
    double sigma0;
    double tmax;
    double xc;

    double Gamma0;
    double x_gamma;
    double d;
    int n;

    int N = 1+L/dx;
    double a=L/10;

};

struct Tridiagonal{
    std::vector<std::complex<double>> diag, e;
};

Init init={
    0.01,    // dx
    0.0001,  // dt
    26.0,    // L ]
    16.0,    // V0 ]
    6.0,    // k0]
    1.2,     // sigma0 ]
    2.8,     // tmax)
    6.5,      // xc
    100.0,     // gamma0
    18,     // x_gamma 
    10.0,     // d
    4        // n
};

void V_init(std::vector<std::complex<double>> &V, std::string pot_type)
{
    if (pot_type=="zero")
    {
        for (auto &v:V) 
            v=0.0;
    }

    else if (pot_type=="imag")
    {
        for (int x_idx=0; x_idx<init.N; x_idx++)
        {
            double x=x_idx*init.dx;

            if (x>init.x_gamma /*|| x<init.L-init.x_gamma*/)
            {
                V[x_idx]=-i*init.Gamma0*pow((x-init.x_gamma)/init.d, init.n);
            }
            else
            {
                V[x_idx]=0.0;
            }
        }
    }

    else if (pot_type=="barrier")
    {
        for (int x_idx=0; x_idx<init.N; x_idx++)
        {
            double x = x_idx*init.dx;
            if (std::abs(init.L/2.0 - x)<init.a/2)
            {
                V[x_idx]=init.V0;
            }
            else if (x>init.x_gamma /*|| x < init.L - init.x_gamma*/)
            {
                V[x_idx]=-i*init.Gamma0*pow((x-init.x_gamma)/init.d, init.n);
            }
            else
            {
                V[x_idx]=0.0;
            }
        }
    }
    
    
}

void Psi_init(std::vector<std::complex<double>> &Psi)
{
    double A = pow(1.0/(2.0*init.sigma0*M_PI), 0.25);
    for (int x_idx=0; x_idx<init.N; x_idx++)
    {
        double x = x_idx*init.dx;
        double x_c = init.xc;
        Psi[x_idx] = A*exp(-(x-x_c)*(x-x_c)/(4.0*init.sigma0*init.sigma0))*exp(i*init.k0*(x-x_c));
    }
}


Tridiagonal Hamiltionian(std::vector<std::complex<double>> &V)
{
    Tridiagonal H;
    H.diag.resize(init.N);
    H.e.resize(init.N);

    double cst=-h*h/(2.0*m);

    for (int I=0; I<init.N; I++)
    {
        H.diag[I]=cst*(-2.0)/(init.dx*init.dx) + V[I];
        H.e[I]=cst/(init.dx*init.dx);
    }

    return H;
}


std::vector<std::complex<double>> b_vec(std::vector<std::complex<double>> &Psi, const Tridiagonal &H)
{
    std::vector<std::complex<double>> b(init.N, 0.0);
    for (int I=1; I<init.N-1; I++)
    {
        b[I]=Psi[I]+0.5*i*init.dt*((-H.diag[I])*Psi[I] - H.e[I]*(Psi[I-1]+Psi[I+1]));
    }

    b[0]=Psi[0] + 0.5*i*init.dt*((-H.diag[0])*Psi[0] - H.e[0]*Psi[1]);
    b[init.N-1]=Psi[init.N-1] + 0.5*i*init.dt*((-H.diag[init.N-1])*Psi[init.N-1] - H.e[init.N-2]*Psi[init.N-2]);
    return b;
}


Tridiagonal A_matrix(const Tridiagonal &H)
{
    Tridiagonal A;
    A.diag.resize(init.N);
    A.e.resize(init.N);

    for (int I=0; I<init.N; I++)
    {
        A.diag[I]=1.0 + 0.5*i*init.dt*H.diag[I];
        A.e[I]=0.5*i*init.dt*H.e[I];
    }

    return A;
}

void save(std::vector<std::complex<double>> &VEC, std::string filename)
{
    std::ofstream file(filename);
    for (int I=0; I<init.N; I++)
    {
        double x=I*init.dx;
        file<<x<<" "<<VEC[I].real()<<" "<<VEC[I].imag()<<" "<<std::abs(VEC[I])<<"\n";
    }
}

std::vector<std::complex<double>> thomas(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b, std::vector<std::complex<double>> c, std::vector<std::complex<double>> d)
{
    /*
        Ax=d => a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i  
    */
    
    int n = init.N;
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
    std::vector<std::complex<double>> Psi(init.N, 0.0), V(init.N, 0.0), b(init.N, 0.0);
    //pot_type= zero, imag, barrier

    Psi_init(Psi);
    V_init(V, "barrier");
    Tridiagonal H = Hamiltionian(V);
    Tridiagonal A = A_matrix(H);
    double t=0;
    int idx=0;
    while (t<init.tmax)
    {
        t+=init.dt;
        idx++;

        b=b_vec(Psi, H);
        Psi=thomas(A.e, A.diag, A.e, b);

        if(idx%50==0)
        {
            std::string filename = "./outdir/"+std::to_string(idx)+".dat";
            save(Psi, filename);
        }

    }


    return 0;
}