#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    mat_ctx_t ctx;

    printf("all... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_long(ctx);
    {
        mpoly_t a;

        printf("\n");
        printf("mpoly_set_str: \"3  (3)[1 2 3]\"\n");
        fflush(stdout);

        mpoly_init(a, 3, ctx);
        mpoly_set_str(a, "3  (3)[1 2 3]", ctx);

        printf("a = "), mpoly_print(a, ctx), printf("\n");
        mpoly_clear(a, ctx);
    }
    mat_ctx_clear(ctx);

    printf("... ");

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
