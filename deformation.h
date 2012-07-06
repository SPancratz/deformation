#ifndef DEFORMATION_H
#define DEFORMATION_H

#include <stdlib.h>
#include <mpir.h>

#include "generics.h"
#include "fmpz.h"
#include "mpoly.h"
#include "mat.h"

#include "gmconnection.h"

typedef struct {
    long N0;
    long N1;
    long N2;
    long N3;
    long K;
    long m;
    long r, s;
} prec_struct;

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

void deformation_precisions(prec_struct *prec, 
                            const fmpz_t p, long a, long n, long d, long degR);

void frob(const mpoly_t P, const fmpz_t t1, 
          const ctx_t ctxFracQt, const fmpz_t p);

void frob_with_precisions(mat_t F, const ctx_t ctxF, 
                          const mpoly_t P, const ctx_t ctxFracQt, 
                          long NWork, long Kfinite, long Kinfinite);

void frob_with_precisions_fmpq(mat_t F, const ctx_t ctxF, 
                               const mpoly_t P, const ctx_t ctxFracQt, 
                               long NWork, long K);

#endif

