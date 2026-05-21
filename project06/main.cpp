#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <complex>


const double hbar=1;
const double m=1;
const double cst=-hbar*hbar/(2*m);
const std::complex<double> I ={0.0, 1.0}; 

struct Init{
    double Lx;
    double Ly;
    double x0;
    double y0;
    int nX;
    int nY;
    double alpha;
    double V0;
    double e;

    double dx=Lx/(nX-1);
    double dy=Ly/(nY-1);

};

Init init={
    20.0, //Lx
    10.0, //Ly
    10.0, //x0
    5.0, //y0
    800, //nX
    400, //nY
    0.0001, //alpha
    0.05, //V0
    1.0 //e
};
//================================= convert to l ==========================
int cTl(int i, int j)
{
    return j+i*init.nY;
}

//=================================== POTENTIAL ===========================
void V_init(std::vector<double> &V)
{
    for (int x_idx=0; x_idx<init.nX; x_idx++)
    {
        for (int y_idx=0; y_idx<init.nY; y_idx++)
        {
            int l = cTl(x_idx,y_idx);
            double x=x_idx*init.dx;
            double y = y_idx*init.dy;
            double D=cosh(init.alpha*(x-init.x0));
            V[l]=(y-init.y0)*(y-init.y0)*init.V0/(D*D);
        }
    }
}
//================================ THOMAS =================================
std::complex<double> thomas(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b, std::vector<std::complex<double>> c, std::vector<std::complex<double>> d)
{
    /*
        Ax=d => a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i  
    */
    
    int n = init.nX;
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

    return x.back();
}
//================================ H_perp solver ===========================
std::vector<double> y_SOLVER()
{
    std::vector<double> epsilon_total(init.nX*init.nY); 

    Eigen::VectorXd H_diag(init.nY);
    Eigen::VectorXd H_e(init.nY-1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> tridiagonalSolver;

    std::vector<double> V(init.nX*init.nY);
    V_init(V);

    double dy2=init.dy*init.dy;
    std::cout<<"DIAGONALIZATION PROCESS\n\n";
    for (int x_idx=0; x_idx<init.nX; x_idx++)
    {
        
        //hamiltionian for given x
        std::cout<<"BUILDING H for x="<<x_idx*init.dx<<"\n";
        for (int y_idx=0; y_idx<init.nY; y_idx++) //co z w.b.
        {
            int l=cTl(x_idx,y_idx);
            H_diag[y_idx]=-2*cst/dy2+V[l];

            if (y_idx==init.nY-1) break;
            H_e[y_idx]=cst/dy2;
        }
        std::cout<<"SOLVING H for x="<<x_idx*init.dx<<"\n";
        tridiagonalSolver.computeFromTridiagonal(H_diag, H_e, Eigen::EigenvaluesOnly);

        for (int y_idx=0; y_idx<init.nY; y_idx++)
        {
            int l = cTl(x_idx,y_idx);
            epsilon_total[l]=tridiagonalSolver.eigenvalues()(y_idx);
        }
    }
    std::cout<<"DIAGONALIZATION PROCESS DONE\n\n";

    return epsilon_total;
}
//========================== QTBM solver =====================================
std::vector<double> QTBM(double E_f, std::vector<double> &epsilon)
{
    std::vector<double> T(init.nY, 0.0);
    std::vector<std::complex<double>> H_diag(init.nX, 0.0);
    std::vector<std::complex<double>> H_e(init.nX, 0.0);
    std::vector<std::complex<double>> b(init.nX, 0.0);
    
    double dx2=init.dx*init.dx;
    
    std::cout<<"QTBM PROCESS START\n\n";
    for (int y_idx=0; y_idx<init.nY; y_idx++)
    {
        std::complex<double> arg = 2.0 * m * (E_f - epsilon[y_idx]) / (hbar * hbar);
        std::complex<double> k = std::sqrt(arg);
        b[0]=-cst*(1.0-exp(2.0*I*k*init.dx))/dx2;

        std::cout<<"BUILDING H FOR eigen group="<<y_idx<<"\n";
        for (int x_idx=0; x_idx<init.nX; x_idx++)
        {
            int l=cTl(x_idx,y_idx);
            //int l =cTl(y_idx,x_idx);
            H_diag[x_idx]=-(2.0*cst/dx2)+epsilon[l]-E_f;
            H_e[x_idx]=cst/dx2;

            if (x_idx==0 || x_idx==init.nX-1)
            {
                std::complex<double> arg = 2.0 * m * (E_f - epsilon[l]) / (hbar * hbar);
                std::complex<double> k = std::sqrt(arg);
                H_diag[x_idx]+=cst*exp(I*k*init.dx)/dx2;
            }
        }
        std::cout<<"SOLVING FOR eigen group ="<<y_idx<<"\n";
        double tau_mod=std::abs(thomas(H_e,H_diag,H_e,b));
        T[y_idx]=tau_mod*tau_mod;
    }
    std::cout<<"QTBM PROCESS DONE\n\n";
    return T;
}
//============================= Conducation ========================
double G(std::vector<double> T)
{
    double G_val=0.0;
    for (int y_idx=0; y_idx<init.nY; y_idx++)
    {
        G_val+=T[y_idx];
    }

    return G_val * (init.e * init.e) / (M_PI * hbar);
}

int main()
{
    std::vector<double> epsilon=y_SOLVER();
    std::vector<double> T(init.nY,0.0);
    std::ofstream file_G("./outdir/G.dat");
    std::ofstream file_V("./outdir/V.dat");
    std::ofstream file_E("./outdir/E.dat");

    std::vector<double> V(init.nX*init.nY);
    V_init(V);
    for (int x_idx=0; x_idx<init.nX; x_idx++)
    {
        for (int y_idx=0; y_idx<init.nY; y_idx++)
        {
            double x = x_idx*init.dx;
            double y = y_idx*init.dy;

            file_V<<x<<" "<<y<<" "<<V[cTl(x_idx,y_idx)]<<"\n";
        }
    }

    double E=0.001;
    double E_max=3;
    double dE=0.001;

    while (E<E_max)
    {
        T=QTBM(E,epsilon);
        file_G<<E<<" "<<G(T)<<"\n";
        
        E+=dE;
    }

    return 0;
}