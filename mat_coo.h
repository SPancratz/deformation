#ifndef MAT_COO_H
#define MAT_COO_H

#include <stdlib.h>
#include <stdio.h>

#include "generics.h"
#include "perm.h"

typedef struct
{
    long m;
    long n;
    long alloc;
    long length;
    char *list;

} mat_coo;

typedef mat_coo mat_coo_t[1];

/* Memory management *********************************************************/

void mat_coo_init(mat_coo_t A, long m, long n, const mat_ctx_t ctx);
void mat_coo_init2(mat_coo_t A, long m, long n, long alloc, const mat_ctx_t ctx);
void mat_coo_realloc(mat_coo_t A, long alloc, const mat_ctx_t ctx);
void mat_coo_clear(mat_coo_t A, int clear, const mat_ctx_t ctx);

void mat_coo_fit_length(mat_coo_t A, long len, const mat_ctx_t ctx);

/* Assignment ****************************************************************/

void mat_coo_set_entry(mat_coo_t A, long i, long j, const void *x, const mat_ctx_t ctx);

void mat_coo_zero(mat_coo_t A, const mat_ctx_t ctx);

/* Randomisation *************************************************************/

void mat_coo_randtest(mat_coo_t A, flint_rand_t state, double d, const mat_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__
int mat_coo_is_zero(const mat_coo_t A, const mat_ctx_t ctx)
{
    return (A->length == 0);
}

/* Input and output **********************************************************/

int mat_coo_debug(const mat_coo_t A, const mat_ctx_t ctx);

int mat_coo_print_dense(const mat_coo_t A, const mat_ctx_t ctx);

#endif

