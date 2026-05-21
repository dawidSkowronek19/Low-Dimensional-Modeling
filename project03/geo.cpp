#include "geo.h"

// ************************ Geometry *******************

Geometry::Geometry(geoInit Parameters, std::string geometry_type):
    m_parameters(Parameters), m_geometry_type(geometry_type)
    {
        std::cout<<"Geometry creation\n";
        int Nx = m_parameters.Nx;
        int Ny = m_parameters.Ny;
        m_rho.resize(Nx*Ny, 0.0);
        m_space.resize(Nx*Ny, 0);
        V_B.resize(Nx*Ny, 0.0);

        //boundary blocking
        for (int x_idx=0; x_idx<Nx; x_idx++)
        {
            int y_idx=0;
            int l=ctL(x_idx, y_idx, Nx);
            m_space[l]=1;
            y_idx=Ny-1;
            l=ctL(x_idx,y_idx, Nx);
            m_space[l]=1;
        }

        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            int x_idx=0;
            int l=ctL(x_idx, y_idx, Nx);
            m_space[l]=1;
            x_idx=Nx-1;
            l=ctL(x_idx,y_idx, Nx);
            m_space[l]=1;
        }
    }

void Geometry::buildGeometry()
{
    int Nx = m_parameters.Nx;
    int Ny = m_parameters.Ny;
    double centerX=m_parameters.Lx/2.0;
    double centerY=m_parameters.Ly/2.0;
    double R = m_parameters.R;

    double dx=m_parameters.dx;
    double dy=m_parameters.dy;
    
    if(m_geometry_type=="circle")
    {
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                int l = ctL(x_idx,y_idx, Nx);

                double x=x_idx*dx, y=y_idx*dy;
                if ((x-centerX)*(x-centerX)+(y-centerY)*(y-centerY)< R*R) {m_rho[l]=m_parameters.rho_cst;}
                else {m_rho[l]=0.0;}
                V_B[l]=0.0;
            }
        }
    }
    
    else if (m_geometry_type=="quantum_wire")
    {
        
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            double y=y_idx*dy;
            for(int x_idx=0; x_idx<Nx; x_idx++)
            {
                int l = ctL(x_idx, y_idx, Nx);
                

                if (y<m_parameters.d_bound)
                {
                    V_B[l] = m_parameters.Vd;
                    m_space[l]=1;
                }

                else if (y>m_parameters.u_bound)
                {
                    V_B[l] = m_parameters.Vu;
                    m_space[l]=1; 
                }

                else
                {
                    m_rho[l]=m_parameters.rho_cst;
                }
            }

            if (y>m_parameters.d_bound && y<m_parameters.u_bound)
            {
                int x_idx_tmp=0;
                int l = ctL(x_idx_tmp, y_idx, Nx);
                m_space[l]=0;

                x_idx_tmp=Nx-1;
                l=ctL(x_idx_tmp, y_idx, Nx);
                m_space[l]=0;
                
            }
            
        }
    }

    else if (m_geometry_type=="quantum_dot")
    {
        for (int y_idx=0; y_idx<Ny; y_idx++)
        {
            double y=y_idx*dy;
            for (int x_idx=0; x_idx<Nx; x_idx++)
            {
                double x= x_idx*dx;
                int l = ctL(x_idx,y_idx, Nx);
                m_rho[l]=m_parameters.rho_cst;
                for (int rod_idx=0; rod_idx<6; rod_idx++)
                {
                    Rod acc_rod=m_parameters.rod[rod_idx];
                    if (std::abs(y-acc_rod.y_c)<acc_rod.height/2.0 && std::abs(x-acc_rod.x_c)<acc_rod.width/2.0)
                    {
                        m_space[l]=1;
                        V_B[l]=m_parameters.Vg;
                        m_rho[l]=0.0;
                        break;
                    }
                }
            }
        }
    }

    
        
}

double Geometry::GetRho(int l) {return m_rho[l];}
double Geometry::GetVB(int l) {return V_B[l];}
int Geometry::GetSpace(int l ) {return m_space[l];}

Geometry::~Geometry(){std::cout<<"Geo destroyed\n";}



//*******************************************************