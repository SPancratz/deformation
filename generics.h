#ifndef GENERICS_H
#define GENERICS_H

#include "flint.h"
#include "long_extras.h"

typedef struct
{
    size_t size;

    void (*init)(void *op);
    void (*clear)(void *op);

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

/* long **********************************************************************/

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

/* mpq_t **********************************************************************/

static void _mpq_init(void *op)
    { mpq_init(op); }
static void _mpq_clear(void *op)
    { mpq_clear(op); }
static void _mpq_set(void *rop, const void *op)
    { mpq_set(rop, op); }
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

static void mat_ctx_init_mpq(mat_ctx_t ctx)
{
    ctx->size              = sizeof(__mpq_struct);

    ctx->init              = &_mpq_init;
    ctx->clear             = &_mpq_clear;
    ctx->set               = &_mpq_set;
    ctx->swap              = &_mpq_swap;
    ctx->zero              = &_mpq_zero;
    ctx->one               = &_mpq_one;
    ctx->randtest          = &_mpq_randtest;
    ctx->randtest_not_zero = &_mpq_randtest_not_zero;
    ctx->equal             = &_mpq_equal;
    ctx->is_zero           = &_mpq_is_zero;
    ctx->is_one            = &_mpq_is_one;
    ctx->add               = &_mpq_add;
    ctx->sub               = &_mpq_sub;
    ctx->mul               = &_mpq_mul;
    ctx->div               = &_mpq_div;
    ctx->print             = &_mpq_print;
}

static void mat_ctx_clear(mat_ctx_t ctx)
    { }

#endif

