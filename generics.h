#ifndef GENERICS_H
#define GENERICS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpir.h>
#include "flint.h"
#include "long_extras.h"
#include "ulong_extras.h"

#include "fmpz_poly_q.h"

typedef struct
{
    size_t size;

    void (*init)(void *op);
    void (*clear)(void *op);

    void (*set)(void *rop, const void *op);
    void (*set_si)(void *rop, long op);
    void (*swap)(void *op1, void *op2);
    void (*zero)(void *rop);
    void (*one)(void *rop);

    void (*randtest)(void *rop, flint_rand_t state);
    void (*randtest_not_zero)(void *rop, flint_rand_t state);

    int (*equal)(const void *op1, const void *op2);
    int (*is_zero)(const void *op);
    int (*is_one)(const void *op);

    void (*neg)(void *rop, const void *op);

    void (*add)(void *rop, const void *op1, const void *op2);
    void (*sub)(void *rop, const void *op1, const void *op2);
    void (*mul)(void *rop, const void *op1, const void *op2);
    void (*div)(void *rop, const void *op1, const void *op2);

    void (*derivative)(void *rop, const void *op);

    int (*print)(const void *op);
    char * (*get_str)(const void *op);
    int (*set_str)(void *rop, const char * str);

} __mat_ctx_struct;

typedef __mat_ctx_struct mat_ctx_t[1];

/* Predefined contexts *******************************************************/

/* long **********************************************************************/

static void ld_init(void *op)
    { *(long *) op = 0; }
static void ld_clear(void *op) { }
static void ld_set(void *rop, const void *op)
    { *(long *) rop = *(const long *) op; }
static void ld_set_si(void *rop, long op)
    { *(long *) rop = op; }
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
static void ld_neg(void *rop, const void *op)
    { *(long *) rop = - (*(long *) op); }
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
static char * ld_get_str(const void *op)
{
    char *s;

    s = malloc(21);
    if (!s)
    {
        printf("ERROR (ld_get_str).\n");
        abort();
    }
    sprintf(s, "%ld", *(long *) op);
    return s;
}
static int ld_set_str(void *rop, const char *str)
{
    *(long *) rop = atoi(str);
    return 1;  /* TODO */
}

static void mat_ctx_init_long(mat_ctx_t ctx)
{
    ctx->size              = sizeof(long);

    ctx->init              = &ld_init;
    ctx->clear             = &ld_clear;
    ctx->set               = &ld_set;
    ctx->set_si            = &ld_set_si;
    ctx->swap              = &ld_swap;
    ctx->zero              = &ld_zero;
    ctx->one               = &ld_one;
    ctx->randtest          = &ld_randtest;
    ctx->randtest_not_zero = &ld_randtest_not_zero;
    ctx->equal             = &ld_equal;
    ctx->is_zero           = &ld_is_zero;
    ctx->is_one            = &ld_is_one;
    ctx->neg               = &ld_neg;
    ctx->add               = &ld_add;
    ctx->sub               = &ld_sub;
    ctx->mul               = &ld_mul;
    ctx->div               = &ld_div;
    ctx->derivative        = NULL;
    ctx->print             = &ld_print;
    ctx->get_str           = &ld_get_str;
    ctx->set_str           = &ld_set_str;
}

/* mpq_t **********************************************************************/

static void _mpq_init(void *op)
    { mpq_init(op); }
static void _mpq_clear(void *op)
    { mpq_clear(op); }
static void _mpq_set(void *rop, const void *op)
    { mpq_set(rop, op); }
static void _mpq_set_si(void *rop, long op)
    { mpq_set_si(rop, op, 1); }
static void _mpq_swap(void *op1, void *op2)
    { mpq_swap(op1, op2); }
static void _mpq_zero(void *rop)
    { mpq_set_ui(rop, 0, 1); }
static void _mpq_one(void *rop)
    { mpq_set_ui(rop, 1, 1); }
static void _mpq_randtest(void *rop, flint_rand_t state)
{
    mpz_rrandomb(mpq_numref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    while (mpz_sgn(mpq_denref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_denref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    mpq_canonicalize(rop);
}
static void _mpq_randtest_not_zero(void *rop, flint_rand_t state)
{
    while (mpz_sgn(mpq_numref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_numref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    while (mpz_sgn(mpq_denref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_denref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    mpq_canonicalize(rop);
}
static int _mpq_equal(const void *op1, const void *op2)
    { return mpq_equal(op1, op2); }
static int _mpq_is_zero(const void *op)
    { return mpq_sgn((__mpq_struct *) op) == 0; }
static int _mpq_is_one(const void *op)
    { return mpq_cmp_ui((__mpq_struct *) op, 1, 1) == 0; }
static void _mpq_neg(void *rop, const void *op)
    { mpq_neg(rop, op); }
static void _mpq_add(void *rop, const void *op1, const void *op2)
    { mpq_add(rop, op1, op2); }
static void _mpq_sub(void *rop, const void *op1, const void *op2)
    { mpq_sub(rop, op1, op2); }
static void _mpq_mul(void *rop, const void *op1, const void *op2)
    { mpq_mul(rop, op1, op2); }
static void _mpq_div(void *rop, const void *op1, const void *op2)
    { mpq_div(rop, op1, op2); }
static int _mpq_print(const void *op)
    { return gmp_printf("%Qd", op); }
static char * _mpq_get_str(const void *op)
    { return mpq_get_str(NULL, 10, op); }
static int _mpq_set_str(void *rop, const char *str)
    { return (mpq_set_str(rop, str, 10) == 0) ? 1 : 0; }

static void mat_ctx_init_mpq(mat_ctx_t ctx)
{
    ctx->size              = sizeof(__mpq_struct);

    ctx->init              = &_mpq_init;
    ctx->clear             = &_mpq_clear;
    ctx->set               = &_mpq_set;
    ctx->set_si            = &_mpq_set_si;
    ctx->swap              = &_mpq_swap;
    ctx->zero              = &_mpq_zero;
    ctx->one               = &_mpq_one;
    ctx->randtest          = &_mpq_randtest;
    ctx->randtest_not_zero = &_mpq_randtest_not_zero;
    ctx->equal             = &_mpq_equal;
    ctx->is_zero           = &_mpq_is_zero;
    ctx->is_one            = &_mpq_is_one;
    ctx->neg               = &_mpq_neg;
    ctx->add               = &_mpq_add;
    ctx->sub               = &_mpq_sub;
    ctx->mul               = &_mpq_mul;
    ctx->div               = &_mpq_div;
    ctx->derivative        = NULL;
    ctx->print             = &_mpq_print;
    ctx->get_str           = &_mpq_get_str;
    ctx->set_str           = &_mpq_set_str;
}

/* Rational functions ********************************************************/

static void _fmpz_poly_q_init(void *rop)
    { fmpz_poly_q_init(rop); }
static void _fmpz_poly_q_clear(void *rop)
    { fmpz_poly_q_clear(rop); }

static void _fmpz_poly_q_set(void *rop, const void *op)
    { fmpz_poly_q_set(rop, op); }
static void _fmpz_poly_q_set_si(void *rop, long op)
    { fmpz_poly_q_set_si(rop, op); }
static void _fmpz_poly_q_swap(void *op1, void *op2)
    { fmpz_poly_q_swap(op1, op2); }
static void _fmpz_poly_q_zero(void *rop)
    { fmpz_poly_q_zero(rop); }
static void _fmpz_poly_q_one(void *rop)
    { fmpz_poly_q_one(rop); }

static void _fmpz_poly_q_randtest(void *rop, flint_rand_t state)
{
    fmpz_poly_q_randtest(rop, state, n_randint(state, 50), 80, 
                                     n_randint(state, 49) + 1, 80);
}
static void _fmpz_poly_q_randtest_not_zero(void *rop, flint_rand_t state)
{
    fmpz_poly_q_randtest_not_zero(rop, state, n_randint(state, 49) + 1, 80, 
                                              n_randint(state, 49) + 1, 80);
}

static int _fmpz_poly_q_equal(const void *op1, const void *op2)
    { return fmpz_poly_q_equal(op1, op2); }
static int _fmpz_poly_q_is_zero(const void *op)
    { return fmpz_poly_q_is_zero(op); }
static int _fmpz_poly_q_is_one(const void *op)
    { return fmpz_poly_q_is_one(op); }

static void _fmpz_poly_q_neg(void *rop, const void *op)
    { fmpz_poly_q_neg(rop, op); }

static void _fmpz_poly_q_add(void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_add(rop, op1, op2); }
static void _fmpz_poly_q_sub(void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_sub(rop, op1, op2); }
static void _fmpz_poly_q_mul(void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_mul(rop, op1, op2); }
static void _fmpz_poly_q_div(void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_div(rop, op1, op2); }

static void _fmpz_poly_q_derivative(void *rop, const void *op)
    { fmpz_poly_q_derivative(rop, op); }

static int _fmpz_poly_q_print(const void *op)
    { return fmpz_poly_q_print(op); }
static char * _fmpz_poly_q_get_str(const void *op)
    { return fmpz_poly_q_get_str(op); }
static int _fmpz_poly_q_set_str(void *rop, const char *str)
    { return fmpz_poly_q_set_str(rop, str); }

static void mat_ctx_init_fmpz_poly_q(mat_ctx_t ctx)
{
    ctx->size              = sizeof(fmpz_poly_q_struct);

    ctx->init              = &_fmpz_poly_q_init;
    ctx->clear             = &_fmpz_poly_q_clear;
    ctx->set               = &_fmpz_poly_q_set;
    ctx->set_si            = &_fmpz_poly_q_set_si;
    ctx->swap              = &_fmpz_poly_q_swap;
    ctx->zero              = &_fmpz_poly_q_zero;
    ctx->one               = &_fmpz_poly_q_one;
    ctx->randtest          = &_fmpz_poly_q_randtest;
    ctx->randtest_not_zero = &_fmpz_poly_q_randtest_not_zero;
    ctx->equal             = &_fmpz_poly_q_equal;
    ctx->is_zero           = &_fmpz_poly_q_is_zero;
    ctx->is_one            = &_fmpz_poly_q_is_one;
    ctx->neg               = &_fmpz_poly_q_neg;
    ctx->add               = &_fmpz_poly_q_add;
    ctx->sub               = &_fmpz_poly_q_sub;
    ctx->mul               = &_fmpz_poly_q_mul;
    ctx->div               = &_fmpz_poly_q_div;
    ctx->derivative        = &_fmpz_poly_q_derivative;
    ctx->print             = &_fmpz_poly_q_print;
    ctx->get_str           = &_fmpz_poly_q_get_str;
    ctx->set_str           = &_fmpz_poly_q_set_str;
}

static void mat_ctx_clear(mat_ctx_t ctx)
    { }

#endif

