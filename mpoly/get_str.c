#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpoly.h"

/**
 * Returns a string representation of the monomial \c x in \f$n+1\f$ variables 
 * as a space-separated list of exponents.
 */
static char * _mpoly_mon_get_str(mon_t x, int n)
{
    char *s, *str;
    exp_t e;
    int i;

    str = malloc(n * 4 + 1);
    if (str == NULL)
    {
        printf("ERROR (_mpoly_mon_get_str).  Memory allocation failed.\n");
        abort();
    }
    
    /* Case i == 0 */
    e = x & MON_BITMASK_BLOCK;
    s = str + sprintf(str, "%lu", e);
    x >>= MON_BITS_PER_EXP;
    
    for (i = 1; i < n; i++)
    {
        e = x & MON_BITMASK_BLOCK;
        s += sprintf(s, " %lu", e);
        x >>= MON_BITS_PER_EXP;
    }
    
    return str;
}

char * mpoly_get_str(const mpoly_t op, const mat_ctx_t ctx)
{
    char *s, *str;
    char **coeffs, **mons;
    int bound, i, num;

    mpoly_iter_t iter;
    mpoly_term t;

    if (mpoly_is_zero(op, ctx))
    {
        str = malloc(2);
        sprintf(str, "%ld" , op->n);
        return str;
    }

    num = RBTREE_SIZE(mpoly, op->dict);

    coeffs = malloc(num * sizeof(char *));
    mons   = malloc(num * sizeof(char *));

    i = 0;
    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        if (ctx->is_one(t->val))
            coeffs[i] = NULL;
        else
            coeffs[i] = ctx->get_str(t->val);
        mons[i] = _mpoly_mon_get_str(t->key, op->n);
        i++;
    }
    mpoly_iter_clear(iter);

    bound = 3;
    for (i = 0; i < num; i++)
    {
        bound += (coeffs[i]) ? strlen(coeffs[i]) : 0;
        bound += strlen(mons[i]);
        bound += 5;  /* For '(', ')', '[', ']', ' ' */
    }
    
    str = malloc(bound);
    s = str + sprintf(str, "%ld", op->n);
    *s++ = ' ';
    *s++ = ' ';

    for (i = 0; i < num; i++)
    {
        /* Add the coefficient */
        if (coeffs[i])
        {
            *s++ = '(';
            s += sprintf(s, "%s", coeffs[i]);
            *s++ = ')';
        }
        
        /* Add the monomial */
        *s++ = '[';
        s += sprintf(s, "%s", mons[i]);
        *s++ = ']';
        
        /* Add the space if necessary, or complete the string */
        *s++ = (i != num - 1) ? ' ' : '\0';
    }
    
    /* Clear up the temporary strings */

    for (i = 0; i < num; i++)
    {
        if (coeffs[i])
            free(coeffs[i]);
        free(mons[i]);
    }
    free(coeffs);
    free(mons);

    return str;
}

