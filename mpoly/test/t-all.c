#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

#include "mpoly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ctx_t ctx;

    printf("all... ");
    fflush(stdout);

    _randinit(state);

    ctx_init_long(ctx);
    {
        mpoly_t a, b;
        char *s;

        mpoly_init(a, 3, ctx);
        mpoly_init(b, 3, ctx);

        mpoly_set_str(a, "3  (3)[1 2 3]", ctx);
        s = mpoly_get_str(a, ctx);
        mpoly_set_str(b, s, ctx);

        result = mpoly_equal(a, b, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), mpoly_print(a, ctx), printf("\n");
            printf("b = "), mpoly_print(b, ctx), printf("\n");
            printf("s = %s\n", s);
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        free(s);
    }
    ctx_clear(ctx);

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
