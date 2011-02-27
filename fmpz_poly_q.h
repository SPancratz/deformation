/******************************************************************************

    Copyright (C) 209, 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_POLY_Q_H
#define FMPZ_POLY_Q_H

#include <stdlib.h>

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

void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop);

/* Memory management *********************************************************/

void fmpz_poly_q_init(fmpz_poly_q_t rop);

void fmpz_poly_q_clear(fmpz_poly_q_t rop);

/* Randomisation *************************************************************/

void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          long len1, mp_bitcnt_t bits1, 
                          long len2, mp_bitcnt_t bits2);

void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   long len1, mp_bitcnt_t bits1, 
                                   long len2, mp_bitcnt_t bits2);

/* Assignment ****************************************************************/

void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

void fmpz_poly_q_set_si(fmpz_poly_q_t rop, long op);

void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2);

static __inline__
void fmpz_poly_q_zero(fmpz_poly_q_t rop)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_si(rop->den, 1);
}

static __inline__ 
void fmpz_poly_q_one(fmpz_poly_q_t rop)
{
    fmpz_poly_set_si(rop->num, 1);
    fmpz_poly_set_si(rop->den, 1);
}

static __inline__ 
void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    fmpz_poly_neg(rop->num, op->num);
    fmpz_poly_set(rop->den, op->den);
}

void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

/* Comparison ****************************************************************/

static __inline__ 
int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)
{
    return fmpz_poly_length(op->num) < 1;
}

static __inline__ 
int fmpz_poly_q_is_one(const fmpz_poly_q_t op)
{
    return fmpz_poly_is_one(op->num) && fmpz_poly_is_one(op->den);
}

static __inline__ 
int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    return fmpz_poly_equal(op1->num, op2->num) && fmpz_poly_equal(op1->den, op2->den);
}

/* Addition and subtraction **************************************************/

void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);
void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Scalar multiplication and division ****************************************/

void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, long x);
void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x);
void fmpz_poly_q_scalar_mul_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x);

void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, long x);
void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x);
void fmpz_poly_q_scalar_div_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x);

/* Multiplication and division ***********************************************/

void fmpz_poly_q_mul_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

void fmpz_poly_q_mul(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

void fmpz_poly_q_div(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Powering ******************************************************************/

void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp);

/* Derivative ****************************************************************/

void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

/* Evaluation ****************************************************************/

int fmpz_poly_q_evaluate(mpq_t rop, const fmpz_poly_q_t f, const mpq_t a);

/* Input and output **********************************************************/

int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s);

char * fmpz_poly_q_get_str(const fmpz_poly_q_t op);
char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x);

int fmpz_poly_q_print(const fmpz_poly_q_t op);
int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x);

#endif

