#ifndef VEC_H
#define VEC_H

#include "generics.h"

/* Memory management *********************************************************/

char * _vec_init(long n, const mat_ctx_t ctx);

void _vec_clear(char *vec, long n, const mat_ctx_t ctx);

/* Randomisation *************************************************************/

void _vec_randtest(char *vec, long n, flint_rand_t state, const mat_ctx_t ctx);

void _vec_randtest_not_zero(char *vec, long n, flint_rand_t state, const mat_ctx_t ctx);

/* Assignment ****************************************************************/

void _vec_set(char *vec1, const char *vec2, long n, const mat_ctx_t ctx);

void _vec_zero(char *vec, long n, const mat_ctx_t ctx);

/* Addition ******************************************************************/

void _vec_neg(char *vec1, const char *vec2, long n, const mat_ctx_t ctx);

void _vec_add(char *res, const char *vec1, const char *vec2, long n, 
              const mat_ctx_t ctx);

void _vec_sub(char *res, const char *vec1, const char *vec2, long n, 
              const mat_ctx_t ctx);

/* Scalar multiplication *****************************************************/

void _vec_scalar_mul(char *res, const char *vec, long n, const char *x, 
                     const mat_ctx_t ctx);

void _vec_scalar_div(char *res, const char *vec, long n, const char *x, 
                     const mat_ctx_t ctx);

/* Comparison ****************************************************************/

int _vec_equal(const char *vec1, const char *vec2, long n, 
                const mat_ctx_t ctx);

int _vec_is_zero(const char *vec, long n, const mat_ctx_t ctx);

/* Input and output **********************************************************/

int _vec_print(const char *vec, long n, const mat_ctx_t ctx);

#endif

