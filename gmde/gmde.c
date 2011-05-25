/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "gmde.h"

/*
    Takes a matrix M over "fmpz_poly_q" and splits it up into C/u of type
    "fmpz_poly".
 */
void gmde_convert_matrix(mat_t C, fmpz_poly_t u, const ctx_t ctxC, 
                         const mat_t M, const ctx_t ctxM)
{
    fmpz_poly_t t;
    long i, j;

    assert(C->m == M->m && C->n == M->n);

    fmpz_poly_init(t);

    fmpz_poly_set_ui(u, 1);
    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_lcm(t, u, fmpz_poly_q_denref((fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)));
            fmpz_poly_swap(u, t);
        }

    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_divides(t, 
                              u, fmpz_poly_q_denref((fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)));
            fmpz_poly_mul((fmpz_poly_struct *) mat_entry(C, i, j, ctxC), 
                          fmpz_poly_q_numref((fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)), t);
        }

    fmpz_poly_clear(t);
}

/*
    Given an $n \times n$ matrix $M$ defined over "fmpz_poly_q" 
    finds the matrices $C$ modulo $t^N$.
 */
void gmde_solve_series(fmpq_mat_struct *C, long N, 
                       const mat_t M, const ctx_t ctxM)
{
    const long n = M->m;

    fmpz_poly_t Mden;
    fmpz_t Mden0;
    long lenR;

    long i, j;

    fmpq_mat_struct *B;
    long lenB;

    /* Initialisation */
    
    fmpz_poly_init(Mden);
    fmpz_init(Mden0);

    /*
       Express M as Mnum / Mden and then re-write Mnum (a matrix of 
       polynomials) as b (a list of matrices over Q)
     */
    {
        ctx_t XXX;
        mat_t Mnum;

        ctx_init_fmpz_poly(XXX);
        mat_init(Mnum, n, n, XXX);

        /* Set Mden, lenR, and Mden0 */
        gmde_convert_matrix(Mnum, Mden, XXX, M, ctxM);
        fmpz_poly_get_coeff_fmpz(Mden0, Mden, 0);
        lenR = fmpz_poly_length(Mden);

        /* Let lenB be the maximum length of a polynomial in Mnum */
        lenB = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                lenB = FLINT_MAX(lenB, fmpz_poly_length((fmpz_poly_struct *) mat_entry(Mnum, i, j, XXX)));
        
        /* Rewrite Mnum (a matrix of polynomials) as b (a list of matrices / Q) */
        B = malloc(lenB * sizeof(fmpq_mat_struct));
        for (i = 0; i < lenB; i++)
            fmpq_mat_init(B + i, n, n);
        
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                const fmpz_poly_struct *poly = 
                    (fmpz_poly_struct *) mat_entry(Mnum, i, j, XXX);
                long k;

                for (k = 0; k < fmpz_poly_length(poly); k++)
                {
                    fmpz_set(fmpq_mat_entry_num(B + k, i, j), 
                             fmpz_poly_get_coeff_ptr(poly, k));
                }
            }

        mat_clear(Mnum, XXX);
        ctx_clear(XXX);
    }

    /* Solve the differential system iteratively */
    {
        fmpq_mat_t mat;
        fmpz_t coeff;

        fmpq_mat_init(mat, n, n);
        fmpz_init(coeff);

        fmpq_mat_one(C + 0);

        for (i = 0; i < N - 1; i++)
        {
            fmpq_mat_zero(C + (i + 1));

            j = (0 < i - lenB + 1) ? (i - lenB + 1) : 0;
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + b[i-j] * C[j]; */
                fmpq_mat_mul(mat, B + (i - j), C + j);
                fmpq_mat_add(C + (i + 1), C + (i + 1), mat);
            }

            j = (0 < i - lenR + 1) ? (i - lenR + 1) : 0;
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + r[i-j+1] * j * C[j]; */
                fmpz_poly_get_coeff_fmpz(coeff, Mden, i - j + 1);
                fmpz_mul_ui(coeff, coeff, j);
                fmpq_mat_scalar_mul_fmpz(mat, C + j, coeff);
                fmpq_mat_add(C + (i + 1), C + (i + 1), mat);
            }

            fmpz_mul_si(coeff, Mden0, -(i + 1));
            fmpq_mat_scalar_div_fmpz(C + (i + 1), C + (i + 1), coeff);
        }

        fmpq_mat_clear(mat);
        fmpz_clear(coeff);
    }

    /* Clean-up */
    
    fmpz_poly_clear(Mden);
    fmpz_clear(Mden0);

    for (i = 0; i < lenB; i++)
        fmpq_mat_clear(B + i);
    free(B);
}

/*
    Converts an array of length $N$ of $n \times n$ matrices over 
    the rationals to an $n \times n$ matrix of rational polynomials.

    The matrix $A$ is expected to be a matrix over objects of type 
    \code{fmpq_poly_t}.
 */
void gmde_solve_convert(mat_t A, const ctx_t ctxA, 
                        const fmpq_mat_struct *C, long N)
{
    long i, j, k;

    assert(N > 0);
    assert(A->m == C->r && A->n == C->c);

    for (i = 0; i < A->m; i++)
        for (j = 0; j < A->n; j++)
        {
            ctxA->zero(ctxA, mat_entry(A, i, j, ctxA));

            for (k = N - 1; k >= 0; k--)
            {
                if (fmpq_sgn(fmpq_mat_entry(C + k, i, j)))
                {
                    fmpq_poly_set_coeff_fmpq(
                        (fmpq_poly_struct *) mat_entry(A, i, j, ctxA), 
                        k, fmpq_mat_entry(C + k, i, j));
                }
            }
        }
}

