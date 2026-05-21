#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "param.h"
#include "convTO_l.h"
#include "geo.h"


class Solver{

    public:
        Solver(geoInit m_parameters, Geometry &builtGEO);
        ~Solver();
        int NP(int i, int j);
        
        void GaussSeidel();
        void SOR();
        void save(std::string filename="def_out.dat");

    private:
        geoInit m_parameters;
        Geometry &m_builtGEO;
        std::vector<double> m_V;
        void error2(int  l, double V_old, double &sum);

};

#endif

