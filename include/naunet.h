#ifndef __NAUNET_H__
#define __NAUNET_H__

#include <cvode/cvode.h>                // prototypes for CVODE fcts., consts.
#include <nvector/nvector_cuda.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#include <sundials/sundials_types.h>    // defs. of realtype, sunindextype
#include <sundials/sundials_math.h>     // contains the macros ABS, SUNSQR, EXP

#include "naunet_userdata.h"
#include "naunet_ode.h"
#include "naunet_macros.h"


class Naunet
{

private:
    // UserData *m_data;
    N_Vector m_y;
    SUNMatrix m_a;

    realtype m_atol;
    realtype m_rtol;
    void *m_cvode_mem;
    SUNLinearSolver m_ls;

    cusparseHandle_t m_cusp_handle;
    cusolverSpHandle_t m_cusol_handle;
    
public:
    Naunet();
    ~Naunet();
    int initSolver();
    int resetSolver(int nsystem);
    int solve(realtype *ab, realtype dt, UserData *data);
};

#endif