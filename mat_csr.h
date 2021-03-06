/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
 */

#ifndef MAT_CSR_H
#define MAT_CSR_H

#include <stdlib.h>
#include <stdio.h>

#include "generics.h"
#include "perm.h"
#include "mat.h"

typedef struct
{
    long m;       /* # of rows    */
    long n;       /* # of columns */

    long alloc;   /* Allocated length of a */
    char *x;      /* Values */
    long *j;      /* Column indices */
    long *p;      /* p[i] is the pointer into x (or j) to the first non-zero element in row i */
    long *lenr;   /* Number of non-zero elements in row i */
} __mat_csr_struct;

typedef __mat_csr_struct mat_csr_t[1];

typedef struct
{
    /* Sparse matrix data;  x is only a reference */
    long m;
    long n;
    char *x;
    long *j;
    long *p;
    long *lenr;

    /* Permutations, blocks */
    long *pi;
    long *qi;
    long *B;
    long nb;

    /* Dense data, LUP decomposition */
    char *entries;
    char **LU;
    long *P;
} __mat_csr_solve_struct;

typedef __mat_csr_solve_struct mat_csr_solve_t[1];

/* Memory management *********************************************************/

void mat_csr_init(mat_csr_t A, long m, long n, const ctx_t ctx);

void mat_csr_init2(mat_csr_t A, long m, long n, long alloc, const ctx_t ctx);

void mat_csr_realloc(mat_csr_t A, long alloc, const ctx_t ctx);

void mat_csr_clear(mat_csr_t A, const ctx_t ctx);

void mat_csr_fit_length(mat_csr_t A, long len, const ctx_t ctx);

void _mat_csr_sort_rows(long m, char *x, long *j, long *p, long *lenr, 
                        const ctx_t ctx);

void mat_csr_sort_rows(mat_csr_t A, const ctx_t ctx);

/* Assignment ****************************************************************/

void mat_csr_set_mat(mat_csr_t A, const mat_t mat, const ctx_t ctx);

void mat_csr_set_array3(mat_csr_t A, char *mem, long len, int copy, const ctx_t ctx);

void mat_csr_zero(mat_csr_t A, const ctx_t ctx);

/* Randomisation *************************************************************/

void mat_csr_randtest(mat_csr_t A, flint_rand_t state, double d, const ctx_t ctx);

/* Comparison ****************************************************************/

int mat_csr_is_zero(const mat_csr_t A, const ctx_t ctx);

/* Permutations **************************************************************/

void _mat_csr_permute_rows(long m, long *p, long *lenr, const long *pi);

void mat_csr_permute_rows(mat_csr_t A, const long *pi, const ctx_t ctx);

void _mat_csr_permute_cols(long m, long n, long *j, long *p, long *lenr, const long *pi);

void mat_csr_permute_cols(mat_csr_t A, const long *pi, const ctx_t ctx);

/* Multiplication ************************************************************/

void mat_csr_mul_vec(char *res, const mat_csr_t mat, const char *vec, const ctx_t ctx);

/* Linear systems ************************************************************/

long _mat_csr_zfdiagonal(long *pi, long n, const long *j, const long *p, 
                           const long *lenr, long *w);

long mat_csr_zfdiagonal(long *pi, const mat_csr_t A);

long _mat_csr_block_triangularise(long *arp, long *b, long n, const long *j, 
                                  const long *p, const long *lenr, long *w);

long mat_csr_block_triangularise(long *pi, long *b, const mat_csr_t A, 
                                                    const ctx_t ctx);

void mat_csr_solve_init(mat_csr_solve_t s, const mat_csr_t mat, 
                        const ctx_t ctx);

void mat_csr_solve_clear(mat_csr_solve_t s, const ctx_t ctx);

void mat_csr_solve(char *x, const mat_csr_solve_t s, const char *b, 
                   const ctx_t ctx);

/* Input and output **********************************************************/

int mat_csr_debug(const mat_csr_t A, const ctx_t ctx);

int _mat_csr_print_dense(long m, long n, const char *x, const long *j, 
                                         const long *p, const long *lenr, 
                                         const ctx_t ctx);

int mat_csr_print_dense(const mat_csr_t A, const ctx_t ctx);

#endif

