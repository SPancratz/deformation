/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef DIAGFROB_H
#define DIAGFROB_H

#include <stdlib.h>
#include <stdio.h>

#include "generics.h"
#include "flint.h"
#include "padic.h"

#include "mat.h"

#include "mon.h"

#define DIAGFROB_MOD(x, m)  (((x) % (m) >= 0) ? (x) % (m) : (x) % (m) + (m))

static __inline__ 
long fdiv_si(long n, long d)
{
    return ((n < 0) && (n % d) ? n / d - 1 : n / d);
}

void diagfrob_falling_fac_mpq(mpq_t rop, long u, long d, long r);

void diagfrob_coefficient_mpq(mpq_t rop, long m, long p);

void diagfrob_coefficient(padic_t rop, long m, const padic_ctx_t ctx);

void diagfrob_alpha_mpq(mpq_t rop, fmpz *a, long n, long d, 
                        const mon_t u, const mon_t v, 
                        const padic_ctx_t ctx);

void diagfrob_alpha(padic_t rop, const fmpz *a, long n, long d, 
                    const mon_t u, const mon_t v, 
                    const padic_ctx_t ctx);

#endif

