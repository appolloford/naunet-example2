#include "naunet.h"

// check_flag function is from the cvDiurnals_ky.c example from the CVODE package.
// Check function return value...
//   opt == 0 means SUNDIALS function allocates memory so check if
//            returned NULL pointer
//   opt == 1 means SUNDIALS function returns a flag so check if
//            flag >= 0
//   opt == 2 means function allocates memory so check if returned
//            NULL pointer
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *)flagvalue;
        if (*errflag < 0)
        {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return 1;
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    return 0;
}

Naunet::Naunet() : m_atol(1e-20),
                   m_rtol(1e-5)
{
};

Naunet::~Naunet()
{
    // N_VDestroy(m_y);
    N_VFreeEmpty(m_y);
    SUNMatDestroy(m_a);
    CVodeFree(&m_cvode_mem);
    SUNLinSolFree(m_ls);
    // delete m_data;

    cusparseDestroy(m_cusp_handle);
    cusolverSpDestroy(m_cusol_handle);
    };

int Naunet::initSolver()
{

    cusparseCreate(&m_cusp_handle);
    cusolverSpCreate(&m_cusol_handle);
    m_y = N_VNew_Cuda(MAXNGROUPS*NSPECIES);
    m_a = SUNMatrix_cuSparse_NewBlockCSR(MAXNGROUPS, NSPECIES, NSPECIES, NNZ, m_cusp_handle);
    m_ls = SUNLinSol_cuSolverSp_batchQR(m_y, m_a, m_cusol_handle);
    // abstol = N_VNew_Cuda(neq);
    SUNMatrix_cuSparse_SetFixedPattern(m_a, 1);
    jacInit(m_a);

    
    m_cvode_mem = CVodeCreate(CV_BDF);

    int flag;
    flag = CVodeInit(m_cvode_mem, fex, 0.0, m_y);
    if (check_flag(&flag, "CVodeInit", 1))
        return 1;
    flag = CVodeSStolerances(m_cvode_mem, m_rtol, m_atol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
        return 1;
    flag = CVodeSetLinearSolver(m_cvode_mem, m_ls, m_a);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
        return 1;
    flag = CVodeSetJacFn(m_cvode_mem, jac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
        return 1;

    // reset the n_vector to empty, maybe not necessary
    // N_VDestroy(m_y);
    // m_y = N_VNewEmpty_Cuda();

    return NAUNET_SUCCESS;
};

// To reset the size of cusparse solver
int Naunet::resetSolver(int nsystem)
{

    N_VDestroy(m_y);
    SUNMatDestroy(m_a);
    SUNLinSolFree(m_ls);
    CVodeFree(&m_cvode_mem);

    m_cvode_mem = CVodeCreate(CV_BDF);
    m_y = N_VNew_Cuda(nsystem*NSPECIES);
    m_a = SUNMatrix_cuSparse_NewBlockCSR(nsystem, NSPECIES, NSPECIES, NNZ, m_cusp_handle);
    m_ls = SUNLinSol_cuSolverSp_batchQR(m_y, m_a, m_cusol_handle);
    SUNMatrix_cuSparse_SetFixedPattern(m_a, 1);
    jacInit(m_a);

    int flag;
    flag = CVodeInit(m_cvode_mem, fex, 0.0, m_y);
    if (check_flag(&flag, "CVodeInit", 1))
        return 1;
    flag = CVodeSStolerances(m_cvode_mem, m_rtol, m_atol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
        return 1;
    flag = CVodeSetLinearSolver(m_cvode_mem, m_ls, m_a);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
        return 1;
    flag = CVodeSetJacFn(m_cvode_mem, jac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
        return 1;

    // reset the n_vector to empty, maybe not necessary
    // N_VDestroy(m_y);
    // m_y = N_VNewEmpty_Cuda();

    return NAUNET_SUCCESS;
};

int Naunet::solve(realtype *ab, realtype dt, UserData *data)
{
    int flag;

    // This way is too slow
    // realtype *ydata = N_VGetArrayPointer(m_y);
    // for (int i=0; i<NSPECIES; i++)
    // {
    //     ydata[i] = ab[i];
    // }
    N_VSetHostArrayPointer_Cuda(ab, m_y);
    N_VCopyToDevice_Cuda(m_y);

    
// #ifdef NAUNET_DEBUG
//     sunindextype lrw, liw;
//     N_VSpace_Cuda(m_y, &lrw, &liw);
//     printf("NVector space: real-%d, int-%d\n", lrw, liw);
// #endif

    flag = CVodeReInit(m_cvode_mem, 0.0, m_y);
    if (check_flag(&flag, "CVodeReInit", 1))
        return 1;
    flag = CVodeSetUserData(m_cvode_mem, data);
    if (check_flag(&flag, "CVodeSetUserData", 1))
        return 1;

    realtype t0 = 0.0;
    flag = CVode(m_cvode_mem, dt, m_y, &t0, CV_NORMAL);

    N_VCopyFromDevice_Cuda(m_y);
    ab = N_VGetHostArrayPointer_Cuda(m_y);
    
    return NAUNET_SUCCESS;
};