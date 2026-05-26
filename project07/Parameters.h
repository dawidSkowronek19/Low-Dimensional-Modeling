#ifndef PARAM_H
#define PARAM_H

struct Parameters{
    double a;
    double V0;
    int N_k;
    int N_x;

    double w=a/2.0;


};


const Parameters init={
    1,//a
    100,//V0

    200, //N_k
    200, //N_x
};

#endif