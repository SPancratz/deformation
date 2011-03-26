/******************************************************************************

    Copyright (C) 2009, 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef MPOLY_H
#define MPOLY_H

#include "generics.h"
#include "rbtree.h"
#include "mon.h"

/*
   We define alternative key words for "asm" and "inline", allowing 
   the code to be compiled with the "-ansi" flag under GCC
 */
#ifndef __GNUC__
    #define __asm__     asm
    #define __inline__  inline
#endif

/* Red-black tree ************************************************************/

RBTREE_PROTOTYPE_H(mpoly, mon_t, void *, static)
RBTREE_PROTOTYPE_C(mpoly, mon_t, void *, static)

/* Definitions ***************************************************************/

typedef struct 
{
    RBTREE_T(mpoly) dict;  /* Dictionary */
    long n;                /* Number of variables minus one */
} __mpoly_struct;

typedef __mpoly_struct mpoly_t[1];

typedef RBTREE_NODE(mpoly) * mpoly_term;

typedef RBTREE_ITER_T(mpoly) mpoly_iter_t;

/* Memory management *********************************************************/

void mpoly_init(mpoly_t rop, long n, const mat_ctx_t ctx);

void mpoly_clear(mpoly_t rop, const mat_ctx_t ctx);

/* Iterator ******************************************************************/

void mpoly_iter_init(mpoly_iter_t iter, const mpoly_t op);

void mpoly_iter_clear(mpoly_iter_t iter);

mpoly_term mpoly_iter_next(mpoly_iter_t iter);

/* Randomisation *************************************************************/

void mpoly_randtest(mpoly_t rop, flint_rand_t state, long d, long N, 
                    const mat_ctx_t ctx);

void mpoly_randtest_not_zero(mpoly_t rop, flint_rand_t state, long d, long N, 
                             const mat_ctx_t ctx);

void mpoly_randhom(mpoly_t rop, flint_rand_t state, long d, long N, 
                   const mat_ctx_t ctx);

/* Assignment and basic manipulation *****************************************/

void mpoly_set(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx);

void mpoly_swap(mpoly_t op1, mpoly_t op2, const mat_ctx_t ctx);

void mpoly_zero(mpoly_t rop, const mat_ctx_t ctx);

void mpoly_neg(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx);

/* Setting and retrieving coefficients ***************************************/

void mpoly_set_coeff(mpoly_t rop, const mon_t m, const void *c, 
                     const mat_ctx_t ctx);

void mpoly_get_coeff(void *rop, const mpoly_t op, const mon_t m, 
                     const mat_ctx_t ctx);

/* Properties ****************************************************************/

long mpoly_degree(const mpoly_t op, long var, const mat_ctx_t ctx);

/* Comparison ****************************************************************/

int mpoly_is_zero(const mpoly_t op, const mat_ctx_t ctx);

int mpoly_is_one(const mpoly_t op, const mat_ctx_t ctx);

int mpoly_equal(const mpoly_t op1, const mpoly_t op2, const mat_ctx_t ctx);

/* Addition/ subtraction *****************************************************/

void _mpoly_add_in_place(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx);

void mpoly_add(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const mat_ctx_t ctx);

void _mpoly_sub_in_place(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx);

void mpoly_sub(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const mat_ctx_t ctx);

void mpoly_add_coeff(mpoly_t rop, const mon_t m, const void *x, 
                     const mat_ctx_t ctx);

void mpoly_sub_coeff(mpoly_t rop, const mon_t m, const void *x, 
                     const mat_ctx_t ctx);

/* Scalar multiplication *****************************************************/

void mpoly_scalar_mul(mpoly_t rop, const mpoly_t op, const void *x, 
                      const mat_ctx_t ctx);

void mpoly_scalar_mul_si(mpoly_t rop, const mpoly_t op, long x, 
                         const mat_ctx_t ctx);

void mpoly_scalar_div(mpoly_t rop, const mpoly_t op, const void *x, 
                      const mat_ctx_t ctx);

void mpoly_scalar_div_si(mpoly_t rop, const mpoly_t op, long x, 
                         const mat_ctx_t ctx);

/* Multiplication ************************************************************/

void mpoly_mul_mon(mpoly_t rop, const mpoly_t op, const mon_t m, 
                   const mat_ctx_t ctx);

void mpoly_mul(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const mat_ctx_t ctx);

void mpoly_addmul(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
                  const mat_ctx_t ctx);

void mpoly_submul(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
                  const mat_ctx_t ctx);

/* Derivative ****************************************************************/

void mpoly_derivative(mpoly_t rop, const mpoly_t op, int var, 
                      const mat_ctx_t ctx);

/* Input and output **********************************************************/

int mpoly_set_str(mpoly_t rop, const char *str, const mat_ctx_t ctx);

char * mpoly_get_str(const mpoly_t op, const mat_ctx_t ctx);

static __inline__
void mpoly_print(const mpoly_t op, const mat_ctx_t ctx)
{
    char *str;

    str = mpoly_get_str(op, ctx);
    printf("%s", str);
    free(str);
}

#endif

