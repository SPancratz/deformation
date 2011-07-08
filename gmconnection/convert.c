#include "gmconnection.h"

void gmc_convert(fmpz_poly_mat_t numM, fmpz_poly_t denM, 
                 const mat_t M, const ctx_t ctx)
{
    long i, j;
    fmpz_poly_t t;

    fmpz_poly_init(t);

    fmpz_poly_set_ui(denM, 1);
    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_lcm(t, denM, fmpz_poly_q_denref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctx)));
            fmpz_poly_swap(denM, t);
        }

    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_divides(t, denM, fmpz_poly_q_denref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctx)));
            fmpz_poly_mul(fmpz_poly_mat_entry(numM, i, j), 
                fmpz_poly_q_numref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctx)), t);
        }

    fmpz_poly_clear(t);
}

