#ifndef __NAUNET_USERDATA_H__
#define __NAUNET_USERDATA_H__

// Struct for holding the nessesary additional variables for the problem.
typedef struct
{
    double nH;
    double Tgas;
    double user_crflux;
    double user_Av;
    double user_GtoDN;
    
} UserData;

#endif