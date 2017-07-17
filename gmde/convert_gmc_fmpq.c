/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "flint/fmpz_poly_mat.h"
#include "gmconnection.h"
#include "gmde.h"

void gmde_convert_gmc_fmpq(fmpq_mat_struct **numM, long *lenM, fmpz_poly_t denM, 
                           const mat_t M, const ctx_t ctxM)
{
    long i, j, k, n = M->m;
    fmpz_poly_mat_t t;

    assert(M->m == M->n);

    fmpz_poly_mat_init(t, n, n);
    gmc_convert(t, denM, M, ctxM);

    *lenM = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            *lenM = FLINT_MAX(*lenM, 
                fmpz_poly_length(fmpz_poly_mat_entry(t, i, j)));
    
    *numM = malloc(*lenM * sizeof(fmpq_mat_struct));

    if (*numM == NULL)
    {
        printf("ERROR (gmde_convert_gmc_fmpq).  Cannot allocate memory.\n");
        abort();
    }

    for (i = 0; i < *lenM; i++)
        fmpq_mat_init(*numM + i, n, n);
    
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            const fmpz_poly_struct *poly = fmpz_poly_mat_entry(t, i, j);

            for (k = 0; k < fmpz_poly_length(poly); k++)
            {
                fmpz_set(fmpq_mat_entry_num(*numM + k, i, j), 
                         fmpz_poly_get_coeff_ptr(poly, k));
            }
        }

    fmpz_poly_mat_clear(t);
}

