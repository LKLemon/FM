#ifndef DY4_MODE_H
#define DY4_MODE_H

#include <algorithm>

typedef struct Parameter {
    int RF;
    int IF;
    int AF;
    int D0;
    int D1;
    int U;
    int SPS;
    int GCD_RDS;
    int U_RDS;
    int D_RDS;

} parameter;

parameter set_parameter(int mode);

#endif
