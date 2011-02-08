#ifndef MAT_DENSE_H
#define MAT_DENSE_H

#include "mat.h"

typedef struct
{
    char *entries;
    long m;
    long n;
    char **rows;
} __mat_dense_struct;

typedef __mat_dense_struct mat_dense_t[1];

/* Memory management *********************************************************/

void mat_dense_init(mat_dense_t mat, long m, long n, const mat_ctx_t);

void mat_dense_set(mat_dense_t mat1, const mat_dense_t mat2, 
                   const mat_ctx_t ctx);

void mat_dense_clear(mat_dense_t mat, const mat_ctx_t ctx);

#define mat_dense_entry(mat, i, j, ctx)  ((mat)->rows[i] + (j) * (ctx)->size)

void mat_dense_zero(mat_dense_t mat, const mat_ctx_t ctx);

void mat_dense_one(mat_dense_t mat, const mat_ctx_t ctx);

/* Randomisation *************************************************************/

void mat_dense_randtest(mat_dense_t mat, flint_rand_t state, 
                        const mat_ctx_t ctx);

void mat_dense_randrank(mat_dense_t mat, flint_rand_t state, long rank, 
                        const mat_ctx_t ctx);

void mat_dense_randops(mat_dense_t mat, flint_rand_t state, long count, 
                       const mat_ctx_t ctx);

/* Comparison ****************************************************************/

int mat_dense_equal(const mat_dense_t mat1, const mat_dense_t mat2, 
                    const mat_ctx_t ctx);

int mat_dense_is_one(const mat_dense_t mat, const mat_ctx_t ctx);

int mat_dense_is_zero(const mat_dense_t mat, const mat_ctx_t ctx);

/* Matrix addition ***********************************************************/

void _mat_dense_add(char **rowsC, char ** const rowsA, char **const rowsB, 
                    long m, long n, const mat_ctx_t ctx);

void mat_dense_add(mat_dense_t C, const mat_dense_t A, const mat_dense_t B, 
                   const mat_ctx_t ctx);

/* Matrix-vector multiplication **********************************************/

void _mat_dense_mul_vec(char *y, char ** const rows, long m, long n, 
                        const char *x, const mat_ctx_t ctx);

void mat_dense_mul_vec(char *y, const mat_dense_t A, const char *x, 
                       const mat_ctx_t ctx);

/* Matrix multiplication *****************************************************/

void _mat_dense_mul_classical(char **rowsC, 
                              char ** const rowsA, char ** const rowsB, 
                              long ell, long m, long n, const mat_ctx_t ctx);

void mat_dense_mul_classical(mat_dense_t C, const mat_dense_t A, 
                             const mat_dense_t B, const mat_ctx_t ctx);

/* Permutations **************************************************************/

void _mat_dense_permute_rows(char **rows, long m, const long *pi, char **w);

void mat_dense_permute_rows(mat_dense_t mat, const long *pi, 
                            const mat_ctx_t ctx);

/* Linear systems ************************************************************/

void _mat_dense_lup_solve(char *x, char ** const rows, long m, long n, 
                          const long *pi, const char *b, const mat_ctx_t ctx);

void mat_dense_lup_solve(char *x, const mat_dense_t mat, const long *pi, 
                         const char *b, const mat_ctx_t ctx);

int 
_mat_dense_lup_decompose(long *pi, char **rows, long m, const mat_ctx_t ctx);

int mat_dense_lup_decompose(mat_dense_t out, long *pi, const mat_dense_t mat, 
                            const mat_ctx_t ctx);

/* Input and output **********************************************************/

int mat_dense_debug(const mat_dense_t mat, const mat_ctx_t ctx);

int mat_dense_print(const mat_dense_t mat, const mat_ctx_t ctx);

#endif

