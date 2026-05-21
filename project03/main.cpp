#include "geo.h"
#include "solver.h"

int main()
{
    
    omp_set_num_threads(20);
    
    
    /*
    //###################################### DISK ###################################
    geoInit init_disk={

        0.3, //R
        1.0, //Lx
        1.0, //Ly
        1.0, //rho_cst
        200, //Nx
        200, //Ny
        1.5, //omega
        1e-4,//delta
        "dir",

        -1.0, //Vu
        -2.0, //Vd
        0.0,  //d_bound
        0.0   //u_bound

    };

    for (int i=1; i<8; i++)
    {
        init_disk.Lx=1.0+i*0.5;
        init_disk.Ly=1.0+i*0.5;

        init_disk.dx=init_disk.Lx/(init_disk.Nx-1);
        init_disk.dy=init_disk.Ly/(init_disk.Ny-1);

        std::cout<<"SOLVING DISK FOR L = "<<init_disk.Lx<<"\n";
        Geometry DISK(init_disk, "circle");
        DISK.buildGeometry();
        Solver potential_disk(init_disk, DISK);
        //potential_disk.GaussSeidel();
        potential_disk.SOR();
        potential_disk.save("outdir/SOR_DISK_iterr"+std::to_string(i));
    }
    
    */
    /*
    //################################## WIRE #################################
    geoInit init_qW={

        0.25, //R
        5.0, //Lx
        2.0, //Ly
        10.0, //rho_cst
        500, //Nx
        200, //Ny
        1.5, //omega
        1e-4,//delta
        "BornX",

        -0.5, //Vu
        -1.0, //Vd
        0.25,  //d_bound
        1.75   //u_bound

    };

    //Geometry qW(init_qW, "quantum_wire");
    //qW.buildGeometry();
    //Solver potential_qW(init_qW, qW);
    //potential_qW.GaussSeidel();
    //potential_qW.SOR();
    //potential_qW.save("./outdir/SOR_zad2");
    
    for (int i=1; i<10; i++)
    {
        init_qW.Vd=-1.0-i;

        std::cout<<"SOLVING QUANTUM WIRE FOR Vd = "<<init_qW.Vd<<"\n";
        Geometry qW(init_qW, "quantum_wire");
        qW.buildGeometry();
        Solver potential_qW(init_qW, qW);
        //potential_qW.GaussSeidel();
        potential_qW.SOR();
        potential_qW.save("outdir/SOR_qW_iterr"+std::to_string(i));
    }
    
    */
   
    //################################## DOT ##########################
    double Lx=5, Ly=5;
    Rod r1={Lx/2.0-1.0, 0.78*Ly, 0.44*Ly, Lx/10.0};
    Rod r2={Lx/2.0, 0.82*Ly, 0.36*Ly, Lx/10.0};
    Rod r3={Lx/2.0+1.0, 0.78*Ly, 0.44*Ly, Lx/10.0};

    Rod r4={Lx/2.0-1.0, 0.22*Ly, 0.44*Ly, Lx/10.0};
    Rod r5={Lx/2.0, 0.18*Ly, 0.36*Ly, Lx/10.0};
    Rod r6={Lx/2.0+1.0, 0.22*Ly, 0.44*Ly, Lx/10.0};


    geoInit init_qD={
        0.25, //R
        Lx, //Lx
        Ly, //Ly
        1.0, //rho_cst
        500, //Nx
        500, //Ny
        1.5, //omega
        1e-4,//delta
        "Neuman",

        -1.0, //Vu
        -1.0, //Vd
        0.25,  //d_bound
        0.75,   //u_bound

        -1.0,    //Vg
        {r1, r2, r3, r4, r5, r6}
    };

    //Geometry qD(init_qD, "quantum_dot");
    //qD.buildGeometry();
    //Solver potential_qD(init_qD, qD);
    //potential_qD.GaussSeidel();
    //potential_qD.SOR();
    //potential_qD.save("./outdir/SOR_zad3");

    for (int i=1; i<8; i++)
    {
        //init_qD.rho_cst=i*0.5;
        init_qD.Vg=-i*0.5;
        std::cout<<"SOLVING QD FOR V_rod = "<<init_qD.Vg<<"\n";
        Geometry qD(init_qD, "quantum_dot");
        qD.buildGeometry();
        Solver potential_qD(init_qD, qD);
        //potential_qD.GaussSeidel();
        potential_qD.SOR();
        potential_qD.save("outdir/SOR_quantumDot_iterr"+std::to_string(i));
    }
    
    
    return 0;
}