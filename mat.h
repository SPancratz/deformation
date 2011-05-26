#ifndef MAT_H
#define MAT_H

#include "generics.h"
#include "perm.h"

typedef struct
{
    char *entries;
    long m;
    long n;
    char **rows;
} __mat_struct;

typedef __mat_struct mat_t[1];

/* Memory management *********************************************************/

void mat_init(mat_t mat, long m, long n, const ctx_t ctx);

void mat_set(mat_t mat1, const mat_t mat2, 
                   const ctx_t ctx);

void mat_clear(mat_t mat, const ctx_t ctx);

#define mat_entry(mat, i, j, ctx)  ((mat)->rows[i] + (j) * (ctx)->size)

void mat_zero(mat_t mat, const ctx_t ctx);

void mat_one(mat_t mat, const ctx_t ctx);

/* Randomisation *************************************************************/

void mat_randtest(mat_t mat, flint_rand_t state, 
                        const ctx_t ctx);

void mat_randrank(mat_t mat, flint_rand_t state, long rank, 
                        const ctx_t ctx);

void mat_randops(mat_t mat, flint_rand_t state, long count, 
                       const ctx_t ctx);

/* Comparison ****************************************************************/

int mat_equal(const mat_t mat1, const mat_t mat2, 
                    const ctx_t ctx);

int mat_is_one(const mat_t mat, const ctx_t ctx);

int mat_is_zero(const mat_t mat, const ctx_t ctx);

/* Matrix addition ***********************************************************/

void _mat_add(char **rowsC, char ** const rowsA, char **const rowsB, 
                    long m, long n, const ctx_t ctx);

void mat_add(mat_t C, const mat_t A, const mat_t B, 
                   const ctx_t ctx);

/* Matrix-vector multiplication **********************************************/

void _mat_mul_vec(char *y, char ** const rows, long m, long n, 
                        const char *x, const ctx_t ctx);

void mat_mul_vec(char *y, const mat_t A, const char *x, 
                       const ctx_t ctx);

/* Matrix multiplication *****************************************************/

void _mat_mul_classical(char **rowsC, 
                        char ** const rowsA, char ** const rowsB, 
                        long ell, long m, long n, const ctx_t ctx);

void mat_mul_classical(mat_t C, const mat_t A, 
                                const mat_t B, const ctx_t ctx);

static __inline__ 
void _mat_mul(char **rowsC, char ** const rowsA, char ** const rowsB, 
                            long ell, long m, long n, const ctx_t ctx)
{
    _mat_mul_classical(rowsC, rowsA, rowsB, ell, m, n, ctx);
}

static __inline__ 
void mat_mul(mat_t C, const mat_t A, const mat_t B, const ctx_t ctx)
{
    mat_mul_classical(C, A, B, ctx);
}

/* Permutations **************************************************************/

void _mat_permute_rows(char **rows, long m, const long *pi, char **w);

void mat_permute_rows(mat_t mat, const long *pi, 
                            const ctx_t ctx);

/* Linear systems ************************************************************/

void _mat_lup_solve(char *x, char ** const rows, long m, long n, 
                          const long *pi, const char *b, const ctx_t ctx);

void mat_lup_solve(char *x, const mat_t mat, const long *pi, 
                         const char *b, const ctx_t ctx);

int 
_mat_lup_decompose(long *pi, char **rows, long m, const ctx_t ctx);

int mat_lup_decompose(mat_t out, long *pi, const mat_t mat, 
                            const ctx_t ctx);

/* Matrix inverse ************************************************************/

long mat_inv(mat_t B, const mat_t A, const ctx_t ctx);

/* Charpoly ******************************************************************/

void mat_revcharpoly(char *poly, mat_t mat, const ctx_t ctx);

/* Input and output **********************************************************/

int mat_debug(const mat_t mat, const ctx_t ctx);

int _mat_print(char ** const rows, long m, long n, const ctx_t ctx);

int mat_print(const mat_t mat, const ctx_t ctx);

#endif

