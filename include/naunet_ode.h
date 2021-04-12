#ifndef __NAUNET_ODE_H__
#define __NAUNET_ODE_H__

#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sundials/sundials_types.h>   // defs. of realtype, sunindextype
#include "naunet_macros.h"
#include "naunet_userdata.h"

#define IJth(A,i,j) SM_ELEMENT_D(A,i,j)

__device__ __host__
int calculate_rates(realtype *k, realtype *y, UserData *user_data);
int fex(realtype t, N_Vector u, N_Vector u_dot, void *user_data);
int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, 
        void *user_data, N_Vector tmp);
int jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix jmatrix, 
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int jacInit(SUNMatrix jmatrix);
__global__ void f_kernel(realtype *y, realtype *ydot, UserData *udata, int nsystem);
__global__ void j_kernel(realtype *y, realtype *jdata, UserData *udata, int nsystem);

#endif