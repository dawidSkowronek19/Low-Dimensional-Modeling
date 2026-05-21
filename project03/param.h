#ifndef PARAM_H
#define PARAM_H

struct Rod{
    double x_c=0, y_c=0, height=0, width; //x, y- rod center

};

struct geoInit{
    double R=0, Lx=0, Ly=0, rho_cst=0;
    int Nx=0, Ny=0;
    double omega=0;
    double delta=0;

    std::string boundaryCT="";

    //quantWire
    double Vu=0;
    double Vd=0;
    double d_bound=0;
    double u_bound=0;

    //quantumDot
    double Vg=0;
    Rod rod[6];
    

    double dx=Lx/(Nx-1);
    double dy=Ly/(Ny-1);
    
};

#endif