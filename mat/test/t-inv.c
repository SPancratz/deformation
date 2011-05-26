#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("inv... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 50; i++)
    {
        long n;
        ctx_t ctx;
        mat_t A, B, C;
        long ansB, ansC;

        n = n_randint(state, 20) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, n, n, ctx);
        mat_init(B, n, n, ctx);
        mat_init(C, n, n, ctx);

        mat_randtest(A, state, ctx);

        mat_set(B, A, ctx);
        ansC = mat_inv(C, B, ctx);
        ansB = mat_inv(B, B, ctx);

        result = ((ansB == ansC) && (!ansB || mat_equal(B, C, ctx)));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("A: \n"), mat_print(A, ctx), printf("\n");
            printf("B: \n"), mat_print(B, ctx), printf("\n");
            printf("C: \n"), mat_print(C, ctx), printf("\n");
            printf("ansB = %ld\n", ansB);
            printf("ansC = %ld\n", ansC);
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        ctx_clear(ctx);
    }

    /* Check A * A^{-1} == A^{-1} * A == Id */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 50; i++)
    {
        long n;
        ctx_t ctx;
        mat_t A, B, C, D, I;
        long ansB;

        n = n_randint(state, 20) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, n, n, ctx);
        mat_init(B, n, n, ctx);
        mat_init(C, n, n, ctx);
        mat_init(D, n, n, ctx);
        mat_init(I, n, n, ctx);

        mat_randtest(A, state, ctx);

        ansB = mat_inv(B, A, ctx);
        if (!ansB)
        {
            mat_mul(C, A, B, ctx);
            mat_mul(D, B, A, ctx);
        }

        result = (ansB || (mat_equal(C, D, ctx) && mat_is_one(C, ctx)));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("A: \n"), mat_print(A, ctx), printf("\n");
            printf("B: \n"), mat_print(B, ctx), printf("\n");
            printf("C: \n"), mat_print(C, ctx), printf("\n");
            printf("D: \n"), mat_print(D, ctx), printf("\n");
            printf("ansB = %ld\n", ansB);
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

