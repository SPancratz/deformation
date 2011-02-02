#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>

#include "deformation.h"
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "long_extras.h"

typedef struct
{
    size_t size;

    void (*init)(void *op); /* Used */
    void (*clear)(void *op); /* Used */

    void (*set)(void *rop, const void *op);
    void (*swap)(void *op1, void *op2);
    void (*zero)(void *rop);
    void (*one)(void *rop);

    void (*randtest)(void *rop, flint_rand_t state);
    void (*randtest_not_zero)(void *rop, flint_rand_t state);

    int (*equal)(const void *op1, const void *op2);
    int (*is_zero)(const void *op);
    int (*is_one)(const void *op);

    void (*add)(void *rop, const void *op1, const void *op2);
    void (*sub)(void *rop, const void *op1, const void *op2);
    void (*mul)(void *rop, const void *op1, const void *op2);
    void (*div)(void *rop, const void *op1, const void *op2);

    int (*print)(const void *op);

} mat_ctx;

typedef mat_ctx mat_ctx_t[1];

/* Predefined contexts *******************************************************/

/* long */

static void ld_init(void *op)
    { *(long *) op = 0; }
static void ld_clear(void *op) { }
static void ld_set(void *rop, const void *op)
    { *(long *) rop = *(const long *) op; }
static void ld_swap(void *op1, void *op2)
{
    long t;

    t             = *(long *) op1;
    *(long *) op1 = *(long *) op2;
    *(long *) op2 = t;
}
static void ld_zero(void *rop)
    { *(long *) rop = 0; }
static void ld_one(void *rop)
    { *(long *) rop = 1; }
static void ld_randtest(void *rop, flint_rand_t state)
    { *(long *) rop = z_randtest(state); }
static void ld_randtest_not_zero(void *rop, flint_rand_t state)
    { *(long *) rop = z_randtest_not_zero(state); }
static int ld_equal(const void *op1, const void *op2)
    { return *(long *) op1 == *(long *) op2; }
static int ld_is_zero(const void *op)
    { return *(long *) op == 0; }
static int ld_is_one(const void *op)
    { return *(long *) op == 1; }
static void ld_add(void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 + *(long *) op2; }
static void ld_sub(void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 - *(long *) op2; }
static void ld_mul(void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 * *(long *) op2; }
static void ld_div(void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 / *(long *) op2; }
static int ld_print(const void *op)
    { return printf("%ld", *(long *) op); }

static void mat_ctx_init_long(mat_ctx_t ctx)
{
    ctx->size              = sizeof(long);

    ctx->init              = &ld_init;
    ctx->clear             = &ld_clear;
    ctx->set               = &ld_set;
    ctx->swap              = &ld_swap;
    ctx->zero              = &ld_zero;
    ctx->one               = &ld_one;
    ctx->randtest          = &ld_randtest;
    ctx->randtest_not_zero = &ld_randtest_not_zero;
    ctx->equal             = &ld_equal;
    ctx->is_zero           = &ld_is_zero;
    ctx->is_one            = &ld_is_one;
    ctx->add               = &ld_add;
    ctx->sub               = &ld_sub;
    ctx->mul               = &ld_mul;
    ctx->div               = &ld_div;
    ctx->print             = &ld_print;
}

static void mat_ctx_clear(mat_ctx_t ctx)
    { }

/* Randomisation *************************************************************/

/*
    Generates a random permutation vector of length $n$.

    This function uses the Knuth shuffle to generate a uniformly 
    random permutation without retries.
 */

static void _long_vec_randperm(long *vec, long n, flint_rand_t state)
{
    long i, j, t;

    for (i = 0; i < n; i++)
        vec[i] = i;

    for (i = n - 1; i > 0; i--)
    {
        j = n_randint(state, i + 1);
        t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }
}

static int _long_vec_print(const long *vec, long len)
{
    long i;

    printf("%ld", len);
    if (len > 0)
    {
        printf(" ");
        for (i = 0; i < len; i++)
            printf(" %ld", vec[i]);
    }

    return 1;
}

#endif

