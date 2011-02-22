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

/* Temporary functions *******************************************************/

static __inline__
int fmpz_poly_is_one(const fmpz_poly_t op)
{
    return op->length == 1 && *(op->coeffs) == 1L;
}

static __inline__
int fmpz_poly_is_unit(const fmpz_poly_t op)
{
    return op->length == 1 && (*(op->coeffs) == 1L || *(op->coeffs) == -1L);
}

/* Accessing numerator and denominator ***************************************/

#define fmpz_poly_q_numref(op)  ((op)->num)
#define fmpz_poly_q_denref(op)  ((op)->den)

#endif

