#ifndef DEFORMATION_H
#define DEFORMATION_H

#include <stdlib.h>
#include <mpir.h>

#include "generics.h"
#include "mat.h"
#include "mpoly.h"

#include "fmpz.h"
#include "fmpz_poly.h"
#include "padic_mat.h"
#include "qadic.h"

#include "gmconnection.h"

typedef struct {
    long N0;
    long N1;
    long N2;
    long N3;
    long N3i;
    long N3w;
    long N3iw;
    long N4;
    long K;
    long m;
    long r, s;
} prec_t;

/*
    Extracts diagonal fibre from the multivariate polynomial $P$, 
    which is expected to be at $t = 0$.
 */
static __inline__ 
void mpoly_diagonal_fibre(fmpz *a, const mpoly_t P, const ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    const long d = mpoly_degree(P, -1, ctx);
    const long n = P->n;

    mpq_t x, y;

    mpq_init(x);
    mpq_init(y);

    mpoly_iter_init(iter, P);
    while ((t = mpoly_iter_next(iter)))
    {
        long diag = 1, e, i, var = -1;

        for (i = 0; diag && (i < n); i++)
        {
            e = mon_get_exp(t->key, i);
            if (e)
            {
                if (var >= 0)
                    diag = 0;
                else
                    var = i;
            }
        }

        if (diag)
        {
            fmpz_poly_q_evaluate(y, t->val, x);

            if (mpq_sgn(y))
            {
                fmpz_set_mpz(a + var, mpq_numref(y));
            }
            else
            {
                printf("ERROR (mpoly_diagonal_fibre).  Zero diagonal term.\n");
                abort();
            }
        }
    }
    mpoly_iter_clear(iter);

    mpq_clear(x);
    mpq_clear(y);
}

void deformation_precisions(prec_t *prec, 
                            const fmpz_t p, long a, long n, long d, long degR);

void deformation_revcharpoly(fmpz_poly_t rop, const fmpz_poly_mat_t op, long v, long n, 
                             long N0, long r, long s, const qadic_ctx_t Qq);

void frob(const mpoly_t P, const ctx_t ctxFracQt, 
          const qadic_t t1, const qadic_ctx_t Qq, 
          prec_t *prec, const prec_t *prec_in,
          int verbose);

#endif

