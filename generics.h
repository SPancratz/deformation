#ifndef GENERICS_H
#define GENERICS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpir.h>

#include "flint.h"
#include "long_extras.h"
#include "ulong_extras.h"

#include "fmpq_poly.h"
#include "fmpz_poly_q.h"
#include "padic.h"
#include "padic_poly.h"

typedef struct __ctx_struct 
{
    size_t size;

    void (*init)(const struct __ctx_struct * ctx, void *op);
    void (*clear)(const struct __ctx_struct * ctx, void *op);

    void (*set)(const struct __ctx_struct * ctx, void *rop, const void *op);
    void (*set_si)(const struct __ctx_struct * ctx, void *rop, long op);
    void (*swap)(const struct __ctx_struct * ctx, void *op1, void *op2);
    void (*zero)(const struct __ctx_struct * ctx, void *rop);
    void (*one)(const struct __ctx_struct * ctx, void *rop);

    void (*randtest)(const struct __ctx_struct * ctx, void *rop, flint_rand_t state);
    void (*randtest_not_zero)(const struct __ctx_struct * ctx, void *rop, flint_rand_t state);

    int (*equal)(const struct __ctx_struct * ctx, const void *op1, const void *op2);
    int (*is_zero)(const struct __ctx_struct * ctx, const void *op);
    int (*is_one)(const struct __ctx_struct * ctx, const void *op);

    void (*neg)(const struct __ctx_struct * ctx, void *rop, const void *op);

    void (*add)(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2);
    void (*sub)(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2);
    void (*mul)(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2);
    void (*div)(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2);

    void (*derivative)(const struct __ctx_struct * ctx, void *rop, const void *op);

    int (*print)(const struct __ctx_struct * ctx, const void *op);
    char * (*get_str)(const struct __ctx_struct * ctx, const void *op);
    int (*set_str)(const struct __ctx_struct * ctx, void *rop, const char * str);

    padic_ctx_struct *pctx;

} __ctx_struct;

typedef __ctx_struct ctx_t[1];

/* Predefined contexts *******************************************************/

/* long **********************************************************************/

static void ld_init(const struct __ctx_struct * ctx, void *op)
    { *(long *) op = 0; }
static void ld_clear(const struct __ctx_struct * ctx, void *op) { }
static void ld_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { *(long *) rop = *(const long *) op; }
static void ld_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { *(long *) rop = op; }
static void ld_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
{
    long t;

    t             = *(long *) op1;
    *(long *) op1 = *(long *) op2;
    *(long *) op2 = t;
}
static void ld_zero(const struct __ctx_struct * ctx, void *rop)
    { *(long *) rop = 0; }
static void ld_one(const struct __ctx_struct * ctx, void *rop)
    { *(long *) rop = 1; }
static void ld_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { *(long *) rop = z_randtest(state); }
static void ld_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { *(long *) rop = z_randtest_not_zero(state); }
static int ld_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return *(long *) op1 == *(long *) op2; }
static int ld_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return *(long *) op == 0; }
static int ld_is_one(const struct __ctx_struct * ctx, const void *op)
    { return *(long *) op == 1; }
static void ld_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { *(long *) rop = - (*(long *) op); }
static void ld_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 + *(long *) op2; }
static void ld_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 - *(long *) op2; }
static void ld_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 * *(long *) op2; }
static void ld_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { *(long *) rop = *(long *) op1 / *(long *) op2; }
static int ld_print(const struct __ctx_struct * ctx, const void *op)
    { return printf("%ld", *(long *) op); }
static char * ld_get_str(const struct __ctx_struct * ctx, const void *op)
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
static int ld_set_str(const struct __ctx_struct * ctx, void *rop, const char *str)
{
    *(long *) rop = atoi(str);
    return 1;  /* TODO */
}

static void ctx_init_long(ctx_t ctx)
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

static void _mpq_init(const struct __ctx_struct * ctx, void *op)
    { mpq_init(op); }
static void _mpq_clear(const struct __ctx_struct * ctx, void *op)
    { mpq_clear(op); }
static void _mpq_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { mpq_set(rop, op); }
static void _mpq_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { mpq_set_si(rop, op, 1); }
static void _mpq_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { mpq_swap(op1, op2); }
static void _mpq_zero(const struct __ctx_struct * ctx, void *rop)
    { mpq_set_ui(rop, 0, 1); }
static void _mpq_one(const struct __ctx_struct * ctx, void *rop)
    { mpq_set_ui(rop, 1, 1); }
static void _mpq_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    mpz_rrandomb(mpq_numref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    while (mpz_sgn(mpq_denref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_denref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    mpq_canonicalize(rop);
}
static void _mpq_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    while (mpz_sgn(mpq_numref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_numref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    while (mpz_sgn(mpq_denref((__mpq_struct *) rop)) == 0)
        mpz_rrandomb(mpq_denref((__mpq_struct *) rop), (__gmp_randstate_struct *) state, FLINT_BITS);
    mpq_canonicalize(rop);
}
static int _mpq_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return mpq_equal(op1, op2); }
static int _mpq_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return mpq_sgn((__mpq_struct *) op) == 0; }
static int _mpq_is_one(const struct __ctx_struct * ctx, const void *op)
    { return mpq_cmp_ui((__mpq_struct *) op, 1, 1) == 0; }
static void _mpq_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { mpq_neg(rop, op); }
static void _mpq_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { mpq_add(rop, op1, op2); }
static void _mpq_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { mpq_sub(rop, op1, op2); }
static void _mpq_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { mpq_mul(rop, op1, op2); }
static void _mpq_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { mpq_div(rop, op1, op2); }
static int _mpq_print(const struct __ctx_struct * ctx, const void *op)
    { return gmp_printf("%Qd", op); }
static char * _mpq_get_str(const struct __ctx_struct * ctx, const void *op)
    { return mpq_get_str(NULL, 10, op); }
static int _mpq_set_str(const struct __ctx_struct * ctx, void *rop, const char *str)
    { return (mpq_set_str(rop, str, 10) == 0) ? 1 : 0; }

static void ctx_init_mpq(ctx_t ctx)
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

/* Polynomials over Z ********************************************************/

static void __fmpz_poly_init(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_init(rop); }
static void __fmpz_poly_clear(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_clear(rop); }

static void __fmpz_poly_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_set(rop, op); }
static void __fmpz_poly_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { fmpz_poly_set_si(rop, op); }
static void __fmpz_poly_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { fmpz_poly_swap(op1, op2); }
static void __fmpz_poly_zero(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_zero(rop); }
static void __fmpz_poly_one(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_one(rop); }

static void __fmpz_poly_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpz_poly_randtest(rop, state, n_randint(state, 50), 80);
}
static void __fmpz_poly_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpz_poly_randtest_not_zero(rop, state, n_randint(state, 49) + 1, 80);
}

static int __fmpz_poly_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return fmpz_poly_equal(op1, op2); }
static int __fmpz_poly_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_is_zero((fmpz_poly_struct *) op); }
static int __fmpz_poly_is_one(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_is_one(op); }

static void __fmpz_poly_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_neg(rop, op); }

static void __fmpz_poly_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_add(rop, op1, op2); }
static void __fmpz_poly_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_sub(rop, op1, op2); }
static void __fmpz_poly_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_mul(rop, op1, op2); }
static void __fmpz_poly_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_div(rop, op1, op2); }

static void __fmpz_poly_derivative(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_derivative(rop, op); }

static int __fmpz_poly_print(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_print(op); }
static char * __fmpz_poly_get_str(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_get_str(op); }
static int __fmpz_poly_set_str(const struct __ctx_struct * ctx, void *rop, const char *str)
    { return fmpz_poly_set_str(rop, str); }

static void ctx_init_fmpz_poly(ctx_t ctx)
{
    ctx->size              = sizeof(fmpz_poly_struct);

    ctx->init              = &__fmpz_poly_init;
    ctx->clear             = &__fmpz_poly_clear;
    ctx->set               = &__fmpz_poly_set;
    ctx->set_si            = &__fmpz_poly_set_si;
    ctx->swap              = &__fmpz_poly_swap;
    ctx->zero              = &__fmpz_poly_zero;
    ctx->one               = &__fmpz_poly_one;
    ctx->randtest          = &__fmpz_poly_randtest;
    ctx->randtest_not_zero = &__fmpz_poly_randtest_not_zero;
    ctx->equal             = &__fmpz_poly_equal;
    ctx->is_zero           = &__fmpz_poly_is_zero;
    ctx->is_one            = &__fmpz_poly_is_one;
    ctx->neg               = &__fmpz_poly_neg;
    ctx->add               = &__fmpz_poly_add;
    ctx->sub               = &__fmpz_poly_sub;
    ctx->mul               = &__fmpz_poly_mul;
    ctx->div               = &__fmpz_poly_div;
    ctx->derivative        = &__fmpz_poly_derivative;
    ctx->print             = &__fmpz_poly_print;
    ctx->get_str           = &__fmpz_poly_get_str;
    ctx->set_str           = &__fmpz_poly_set_str;
}

/* Polynomials over Q ********************************************************/

static void __fmpq_poly_init(const struct __ctx_struct * ctx, void *rop)
    { fmpq_poly_init(rop); }
static void __fmpq_poly_clear(const struct __ctx_struct * ctx, void *rop)
    { fmpq_poly_clear(rop); }

static void __fmpq_poly_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpq_poly_set(rop, op); }
static void __fmpq_poly_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { fmpq_poly_set_si(rop, op); }
static void __fmpq_poly_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { fmpq_poly_swap(op1, op2); }
static void __fmpq_poly_zero(const struct __ctx_struct * ctx, void *rop)
    { fmpq_poly_zero(rop); }
static void __fmpq_poly_one(const struct __ctx_struct * ctx, void *rop)
    { fmpq_poly_one(rop); }

static void __fmpq_poly_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpq_poly_randtest(rop, state, n_randint(state, 50), 80);
}
static void __fmpq_poly_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpq_poly_randtest_not_zero(rop, state, n_randint(state, 49) + 1, 80);
}

static int __fmpq_poly_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return fmpq_poly_equal(op1, op2); }
static int __fmpq_poly_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return fmpq_poly_is_zero(op); }
static int __fmpq_poly_is_one(const struct __ctx_struct * ctx, const void *op)
    { return fmpq_poly_is_one(op); }

static void __fmpq_poly_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpq_poly_neg(rop, op); }

static void __fmpq_poly_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpq_poly_add(rop, op1, op2); }
static void __fmpq_poly_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpq_poly_sub(rop, op1, op2); }
static void __fmpq_poly_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpq_poly_mul(rop, op1, op2); }
static void __fmpq_poly_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpq_poly_div(rop, op1, op2); }

static void __fmpq_poly_derivative(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpq_poly_derivative(rop, op); }

static int __fmpq_poly_print(const struct __ctx_struct * ctx, const void *op)
    { return fmpq_poly_print(op); }
static char * __fmpq_poly_get_str(const struct __ctx_struct * ctx, const void *op)
    { return fmpq_poly_get_str(op); }
static int __fmpq_poly_set_str(const struct __ctx_struct * ctx, void *rop, const char *str)
    { return fmpq_poly_set_str(rop, str); }

static void ctx_init_fmpq_poly(ctx_t ctx)
{
    ctx->size              = sizeof(fmpq_poly_struct);

    ctx->init              = &__fmpq_poly_init;
    ctx->clear             = &__fmpq_poly_clear;
    ctx->set               = &__fmpq_poly_set;
    ctx->set_si            = &__fmpq_poly_set_si;
    ctx->swap              = &__fmpq_poly_swap;
    ctx->zero              = &__fmpq_poly_zero;
    ctx->one               = &__fmpq_poly_one;
    ctx->randtest          = &__fmpq_poly_randtest;
    ctx->randtest_not_zero = &__fmpq_poly_randtest_not_zero;
    ctx->equal             = &__fmpq_poly_equal;
    ctx->is_zero           = &__fmpq_poly_is_zero;
    ctx->is_one            = &__fmpq_poly_is_one;
    ctx->neg               = &__fmpq_poly_neg;
    ctx->add               = &__fmpq_poly_add;
    ctx->sub               = &__fmpq_poly_sub;
    ctx->mul               = &__fmpq_poly_mul;
    ctx->div               = &__fmpq_poly_div;
    ctx->derivative        = &__fmpq_poly_derivative;
    ctx->print             = &__fmpq_poly_print;
    ctx->get_str           = &__fmpq_poly_get_str;
    ctx->set_str           = &__fmpq_poly_set_str;
}

/* Rational functions ********************************************************/

static void _fmpz_poly_q_init(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_q_init(rop); }
static void _fmpz_poly_q_clear(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_q_clear(rop); }

static void _fmpz_poly_q_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_q_set(rop, op); }
static void _fmpz_poly_q_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { fmpz_poly_q_set_si(rop, op); }
static void _fmpz_poly_q_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { fmpz_poly_q_swap(op1, op2); }
static void _fmpz_poly_q_zero(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_q_zero(rop); }
static void _fmpz_poly_q_one(const struct __ctx_struct * ctx, void *rop)
    { fmpz_poly_q_one(rop); }

static void _fmpz_poly_q_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpz_poly_q_randtest(rop, state, n_randint(state, 50), 80, 
                                     n_randint(state, 49) + 1, 80);
}
static void _fmpz_poly_q_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
{
    fmpz_poly_q_randtest_not_zero(rop, state, n_randint(state, 49) + 1, 80, 
                                              n_randint(state, 49) + 1, 80);
}

static int _fmpz_poly_q_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return fmpz_poly_q_equal(op1, op2); }
static int _fmpz_poly_q_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_q_is_zero(op); }
static int _fmpz_poly_q_is_one(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_q_is_one(op); }

static void _fmpz_poly_q_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_q_neg(rop, op); }

static void _fmpz_poly_q_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_add(rop, op1, op2); }
static void _fmpz_poly_q_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_sub(rop, op1, op2); }
static void _fmpz_poly_q_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_mul(rop, op1, op2); }
static void _fmpz_poly_q_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { fmpz_poly_q_div(rop, op1, op2); }

static void _fmpz_poly_q_derivative(const struct __ctx_struct * ctx, void *rop, const void *op)
    { fmpz_poly_q_derivative(rop, op); }

static int _fmpz_poly_q_print(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_q_print(op); }
static char * _fmpz_poly_q_get_str(const struct __ctx_struct * ctx, const void *op)
    { return fmpz_poly_q_get_str(op); }
static int _fmpz_poly_q_set_str(const struct __ctx_struct * ctx, void *rop, const char *str)
    { return fmpz_poly_q_set_str(rop, str); }

static void ctx_init_fmpz_poly_q(ctx_t ctx)
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

/* padic_t *******************************************************************/

static void __padic_init(const struct __ctx_struct * ctx, void *op)
    { padic_init(op, ctx->pctx); }
static void __padic_clear(const struct __ctx_struct * ctx, void *op)
    { padic_clear(op, ctx->pctx); }
static void __padic_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { padic_set(rop, op, ctx->pctx); }
static void __padic_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { padic_set_si(rop, op, ctx->pctx); }
static void __padic_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { padic_swap(op1, op2, ctx->pctx); }
static void __padic_zero(const struct __ctx_struct * ctx, void *rop)
    { padic_zero(rop, ctx->pctx); }
static void __padic_one(const struct __ctx_struct * ctx, void *rop)
    { padic_one(rop, ctx->pctx); }
static void __padic_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { padic_randtest(rop, state, ctx->pctx); }
static void __padic_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { padic_randtest_not_zero(rop, state, ctx->pctx); }
static int __padic_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return padic_equal(op1, op2, ctx->pctx); }
static int __padic_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return padic_is_zero(op, ctx->pctx); }
static int __padic_is_one(const struct __ctx_struct * ctx, const void *op)
    { return padic_is_one(op, ctx->pctx); }
static void __padic_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { padic_neg(rop, op, ctx->pctx); }
static void __padic_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_add(rop, op1, op2, ctx->pctx); }
static void __padic_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_sub(rop, op1, op2, ctx->pctx); }
static void __padic_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_mul(rop, op1, op2, ctx->pctx); }
static void __padic_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_div(rop, op1, op2, ctx->pctx); }
static int ___padic_print(const struct __ctx_struct * ctx, const void *op)
    { return padic_print(op, ctx->pctx); }
static char * __padic_get_str(const struct __ctx_struct * ctx, const void *op)
    { return padic_get_str(NULL, op, ctx->pctx); }

static void ctx_init_padic(ctx_t ctx, const padic_ctx_struct * pctx)
{
    ctx->size              = 2 * sizeof(long);

    ctx->init              = &__padic_init;
    ctx->clear             = &__padic_clear;
    ctx->set               = &__padic_set;
    ctx->set_si            = &__padic_set_si;
    ctx->swap              = &__padic_swap;
    ctx->zero              = &__padic_zero;
    ctx->one               = &__padic_one;
    ctx->randtest          = &__padic_randtest;
    ctx->randtest_not_zero = &__padic_randtest_not_zero;
    ctx->equal             = &__padic_equal;
    ctx->is_zero           = &__padic_is_zero;
    ctx->is_one            = &__padic_is_one;
    ctx->neg               = &__padic_neg;
    ctx->add               = &__padic_add;
    ctx->sub               = &__padic_sub;
    ctx->mul               = &__padic_mul;
    ctx->div               = &__padic_div;
    ctx->derivative        = NULL;
    ctx->print             = &___padic_print;
    ctx->get_str           = &__padic_get_str;
    ctx->set_str           = NULL;

    ctx->pctx              = (padic_ctx_struct *) pctx;
}

/* padic_poly_t **************************************************************/

static void __padic_poly_init(const struct __ctx_struct * ctx, void *op)
    { _padic_poly_init(op); }
static void __padic_poly_clear(const struct __ctx_struct * ctx, void *op)
    { _padic_poly_clear(op); }
static void __padic_poly_set(const struct __ctx_struct * ctx, void *rop, const void *op)
    { padic_poly_set(rop, op, ctx->pctx); }
static void __padic_poly_set_si(const struct __ctx_struct * ctx, void *rop, long op)
    { padic_poly_set_si(rop, op, ctx->pctx); }
static void __padic_poly_swap(const struct __ctx_struct * ctx, void *op1, void *op2)
    { _padic_poly_swap(op1, op2); }
static void __padic_poly_zero(const struct __ctx_struct * ctx, void *rop)
    { padic_poly_zero(rop, ctx->pctx); }
static void __padic_poly_one(const struct __ctx_struct * ctx, void *rop)
    { padic_poly_one(rop, ctx->pctx); }
static void __padic_poly_randtest(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { padic_poly_randtest(rop, state, n_randint(state, 50), ctx->pctx); }
static void __padic_poly_randtest_not_zero(const struct __ctx_struct * ctx, void *rop, flint_rand_t state)
    { padic_poly_randtest_not_zero(rop, state, n_randint(state, 49) + 1, ctx->pctx); }
static int __padic_poly_equal(const struct __ctx_struct * ctx, const void *op1, const void *op2)
    { return padic_poly_equal(op1, op2, ctx->pctx); }
static int __padic_poly_is_zero(const struct __ctx_struct * ctx, const void *op)
    { return padic_poly_is_zero(op, ctx->pctx); }
static int __padic_poly_is_one(const struct __ctx_struct * ctx, const void *op)
    { return padic_poly_is_one(op, ctx->pctx); }
static void __padic_poly_neg(const struct __ctx_struct * ctx, void *rop, const void *op)
    { padic_poly_neg(rop, op, ctx->pctx); }
static void __padic_poly_add(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_poly_add(rop, op1, op2, ctx->pctx); }
static void __padic_poly_sub(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_poly_sub(rop, op1, op2, ctx->pctx); }
static void __padic_poly_mul(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
    { padic_poly_mul(rop, op1, op2, ctx->pctx); }
static void __padic_poly_div(const struct __ctx_struct * ctx, void *rop, const void *op1, const void *op2)
{
    printf("ERROR (generics.h).  padic_poly_div not implemented.\n\n");
    abort();
}
static void __padic_poly_derivative(const struct __ctx_struct * ctx, void *rop, const void *op)
    { padic_poly_derivative(rop, op, ctx->pctx); }
static int __padic_poly_print(const struct __ctx_struct * ctx, const void *op)
    { return padic_poly_print(op, ctx->pctx); }
static char * __padic_poly_get_str(const struct __ctx_struct * ctx, const void *op)
    { return padic_poly_get_str(op, ctx->pctx); }

static void ctx_init_padic_poly(ctx_t ctx, const padic_ctx_struct * pctx)
{
    ctx->size              = sizeof(padic_poly_struct);

    ctx->init              = &__padic_poly_init;
    ctx->clear             = &__padic_poly_clear;
    ctx->set               = &__padic_poly_set;
    ctx->set_si            = &__padic_poly_set_si;
    ctx->swap              = &__padic_poly_swap;
    ctx->zero              = &__padic_poly_zero;
    ctx->one               = &__padic_poly_one;
    ctx->randtest          = &__padic_poly_randtest;
    ctx->randtest_not_zero = &__padic_poly_randtest_not_zero;
    ctx->equal             = &__padic_poly_equal;
    ctx->is_zero           = &__padic_poly_is_zero;
    ctx->is_one            = &__padic_poly_is_one;
    ctx->neg               = &__padic_poly_neg;
    ctx->add               = &__padic_poly_add;
    ctx->sub               = &__padic_poly_sub;
    ctx->mul               = &__padic_poly_mul;
    ctx->div               = &__padic_poly_div;
    ctx->derivative        = &__padic_poly_derivative;
    ctx->print             = &__padic_poly_print;
    ctx->get_str           = &__padic_poly_get_str;
    ctx->set_str           = NULL;

    ctx->pctx              = (padic_ctx_struct *) pctx;
}


static void ctx_clear(ctx_t ctx)
    { }

#endif

