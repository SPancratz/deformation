/******************************************************************************

    Copyright (C) 209, 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_POLY_Q_H
#define FMPZ_POLY_Q_H

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

typedef struct
{
    fmpz_poly_struct *num;
    fmpz_poly_struct *den;
} fmpz_poly_q_struct;

typedef fmpz_poly_q_struct fmpz_poly_q_t[1];

#endif

