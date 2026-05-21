#include "solver.h"

// ******************** Constructor/Destructor *************************

Solver::Solver(geoInit init, Geometry &builtGEO): m_parameters(init), m_builtGEO(builtGEO)
{
    std::cout<<"geometry accepted\n";
    m_V.resize(m_parameters.Nx*m_parameters.Ny, 0.0);

    for (int l=0; l<m_parameters.Nx*m_parameters.Ny; l++)
    {
        m_V[l]=m_builtGEO.GetVB(l);
    }
}

Solver::~Solver() {std::cout<<"Solver done\n";}
// **********************************************************************

// ****************************** G-S ***********************************

int Solver::NP(int i, int j){return (i+j)%2==0? 1 : 0;} //return 1- jasne, return 0- ciemne

void Solver::error2(int l ,double V_old, double &sum)
{
    sum+=(m_V[l]-V_old)*(m_V[l]-V_old);
}



void Solver::GaussSeidel()
{
    int Nx = m_parameters.Nx;
    int Ny = m_parameters.Ny;
    double dx=m_parameters.dx;
    double dy=m_parameters.dy;
    
    double errorVal=100;
    std::string bCT = m_parameters.boundaryCT;
    int iterr=0;

    while(errorVal>m_parameters.delta)
    {
        iterr++;
        double error2Val=0.0;
        //std::cout<<"iterr "<<iterr<<"\n";
        #pragma omp parallel for reduction(+:error2Val) collapse(2)
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                double V_old;
                int l =ctL(x_idx, y_idx, Nx);
                if (m_builtGEO.GetSpace(l)) //electrods and dirichlet (not updating nods on boundary and electrods)
                {
                    continue;
                }

                if (NP(x_idx,y_idx)) //day/night nods
                {
                    continue;
                }
                V_old=m_V[l];

                //Born boundary conditions for x only
                if (bCT=="BornX"&&x_idx==0)
                {
                    //std::cout<<"Born at 0\n";
                    int l_born = ctL(Nx-1, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l_born])/(2.0*dx*dx+2.0*dy*dy);
                    error2(l, V_old, error2Val);
                    continue;
                }
                else if (bCT=="BornX"&&x_idx==Nx-1)
                {
                    //std::cout<<"Born at Nx-1\n";
                    int l_born = ctL(0, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l_born]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                    error2(l, V_old, error2Val);
                    continue;
                }
                
                //std::cout<<"solver inside\n";
                m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                error2(l, V_old, error2Val);
            }
        }


        #pragma omp parallel for reduction(+:error2Val) collapse(2)
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                double V_old;
                int l =ctL(x_idx, y_idx, Nx);
                V_old=m_V[l];
                if (m_builtGEO.GetSpace(l)) //electrods and dirichlet (not updating nods on boundary and electrods)
                {
                    continue;
                }

                if (!NP(x_idx,y_idx)) //day/night nods
                {
                    continue;
                }

                //Born boundary conditions for x only
                if (bCT=="BornX"&&x_idx==0)
                {
                    //std::cout<<"Born at 0\n";
                    int l_born = ctL(Nx-1, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l_born])/(2.0*dx*dx+2.0*dy*dy);
                    error2(l, V_old, error2Val);
                    continue;
                }
                else if (bCT=="BornX"&&x_idx==Nx-1)
                {
                    //std::cout<<"Born at Nx-1\n";
                    int l_born = ctL(0, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l_born]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                    error2(l, V_old, error2Val);
                    continue;
                }
                //std::cout<<"solver inside\n";
                m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                error2(l, V_old, error2Val);

            }
        }

        if (bCT=="Neuman")
        {
            for (int y_idx=0; y_idx<Ny-1; y_idx++)
            {
                double V_old;
                int x_idx=0;
                int l=ctL(x_idx,y_idx, Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l+1];
                error2(l, V_old, error2Val);


                x_idx=Nx-1;
                l=ctL(x_idx,y_idx,Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l-1];
                error2(l, V_old, error2Val);
            }

            for (int x_idx=0; x_idx<Nx-1; x_idx++)
            {
                double V_old;
                int y_idx=0;
                int l=ctL(x_idx,y_idx, Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l+Nx];
                error2(l, V_old, error2Val);


                y_idx=Ny-1;
                l=ctL(x_idx,y_idx,Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l-Nx];
                error2(l, V_old, error2Val);
            }
        }
        
        errorVal=sqrt(error2Val);
        
    }
    std::cout<<"iterr "<<iterr<<"\n";  

}

// **********************************************************************


// ******************************** SOR **********************************

void Solver::SOR()
{
    int Nx = m_parameters.Nx;
    int Ny = m_parameters.Ny;
    double dx=m_parameters.dx;
    double dy=m_parameters.dy;
    
    double errorVal=100;
    int iterr=0;
    std::string bCT = m_parameters.boundaryCT;

     while(errorVal>m_parameters.delta)
    {
        iterr++;
        double error2Val=0.0;
        //std::cout<<"iterr "<<iterr<<"\n";

        #pragma omp parallel for reduction(+:error2Val) collapse(2)
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                double V_old;
                int l =ctL(x_idx, y_idx, Nx);
                if (m_builtGEO.GetSpace(l)) //electrods and dirichlet (not updating nods on boundary and electrods)
                {
                    continue;
                }

                if (NP(x_idx,y_idx)) //day/night nods
                {
                    continue;
                }
                V_old=m_V[l];

                //Born boundary conditions for x only
                if (bCT=="BornX"&&x_idx==0)
                {
                    //std::cout<<"Born at 0\n";
                    int l_born = ctL(Nx-1, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l_born])/(2.0*dx*dx+2.0*dy*dy);
                    m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                    error2(l, V_old, error2Val);
                    continue;
                }
                else if (bCT=="BornX"&&x_idx==Nx-1)
                {
                    //std::cout<<"Born at Nx-1\n";
                    int l_born = ctL(0, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l_born]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                    m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                    error2(l, V_old, error2Val);
                    continue;
                }
                
                //std::cout<<"solver inside\n";
                m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);
            }
        }

        #pragma omp parallel for reduction(+:error2Val) collapse(2)
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                double V_old;
                int l =ctL(x_idx, y_idx, Nx);
                V_old=m_V[l];
                if (m_builtGEO.GetSpace(l)) //electrods and dirichlet (not updating nods on boundary and electrods)
                {
                    continue;
                }

                if (!NP(x_idx,y_idx)) //day/night nods
                {
                    continue;
                }

                //Born boundary conditions for x only
                if (bCT=="BornX"&&x_idx==0)
                {
                    //std::cout<<"Born at 0\n";
                    int l_born = ctL(Nx-1, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l_born])/(2.0*dx*dx+2.0*dy*dy);
                    m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                    error2(l, V_old, error2Val);
                    continue;
                }
                else if (bCT=="BornX"&&x_idx==Nx-1)
                {
                    //std::cout<<"Born at Nx-1\n";
                    int l_born = ctL(0, y_idx, Nx);
                    m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l_born]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                    m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                    error2(l, V_old, error2Val);
                    continue;
                }
                //std::cout<<"solver inside\n";
                m_V[l] = (dx*dx*dy*dy*m_builtGEO.GetRho(l)+dx*dx*m_V[l+Nx]+dx*dx*m_V[l-Nx]+dy*dy*m_V[l+1]+dy*dy*m_V[l-1])/(2.0*dx*dx+2.0*dy*dy);
                m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);

            }
        }

        if (bCT=="Neuman")
        {
            for (int y_idx=0; y_idx<Ny-1; y_idx++)
            {
                double V_old;
                int x_idx=0;
                int l=ctL(x_idx,y_idx, Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l+1];
                //m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);


                x_idx=Nx-1;
                l=ctL(x_idx,y_idx,Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l-1];
                //m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);
            }

            for (int x_idx=0; x_idx<Nx-1; x_idx++)
            {
                double V_old;
                int y_idx=0;
                int l=ctL(x_idx,y_idx, Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l+Nx];
                //m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);


                y_idx=Ny-1;
                l=ctL(x_idx,y_idx,Nx);
                V_old=m_V[l];
                m_V[l]=m_V[l-Nx];
                //m_V[l]=m_parameters.omega*m_V[l]+(1.0-m_parameters.omega)*V_old;
                error2(l, V_old, error2Val);
            }
        }
        
        errorVal=sqrt(error2Val);
        //std::cout<<errorVal<<"\n";
    } 
    std::cout<<"iterr "<<iterr<<"\n";
}

// *********************************************************************


// ******************************* SAVE *********************************

void Solver::save(std::string filename)
{
    std::ofstream file(filename);
    for (int y_idx=0; y_idx<m_parameters.Ny; y_idx++)
    {
        for (int x_idx=0; x_idx<m_parameters.Nx; x_idx++)
        {
            file<<x_idx*m_parameters.dx<<" "<<y_idx*m_parameters.dy<<" "<<m_V[ctL(x_idx, y_idx, m_parameters.Nx)]<<"\n";
        }
    }
}

// **********************************************************************