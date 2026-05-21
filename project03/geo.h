#ifndef GEO_H
#define GEO_H

#include <vector>
#include <string>
#include <iostream>
#include "param.h"
#include "convTO_l.h"


class Geometry{
    public:
        Geometry(geoInit parameters, std::string geometry_type);
        ~Geometry();
        void buildGeometry();

        //getters
        double GetRho(int l);
        int GetSpace(int l);
        double GetVB(int l);

    private:
        geoInit m_parameters;
        std::string m_geometry_type;
        std::vector<double> m_rho;
        std::vector<double> V_B;
        std::vector<int> m_space;
};

#endif