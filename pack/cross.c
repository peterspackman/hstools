#include <Python.h>
#include <stdio.h>
#include <math.h>
#include "cio.h"

// COMPUTE CROSS PRODUCT IN C, no error checking

void cross3D(double v1[], double v2[], double r[])
{
    r[0] =   ( (v1[1] * v2[2]) - (v1[2] * v2[1]) );
    r[1] = - ( (v1[0] * v2[2]) - (v1[2] * v2[0]) );
    r[2] =   ( (v1[0] * v2[1]) - (v1[1] * v2[0]) );
}

// COMPUTE NORMAL IN C, no error checking
double normal3D(double v[])
{
    return sqrt( pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2) );
}
