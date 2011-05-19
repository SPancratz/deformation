#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpoly.h"

/*
    Returns the specified substring.

    Given a null-terminated string \code{str} and indices $i$ and $j$, 
    returns the null-terminated substring of \code{str} in the range 
    $[i, j)$.
 */
static char * mpoly_str_substr(const char * str, size_t i, size_t j)
{
    char * substr;
    size_t kstr, ksubstr;
    const size_t len = strlen(str);
    
    if (i < 0)
        i = 0;
    if (j > len)
        j = len;
    
    if (i >= j)
    {
        substr = malloc(1);
        substr[0] = '\0';
        return substr;
    }
    
    substr = malloc(j - i + 1);
    ksubstr = 0;
    for (kstr = i; kstr < j; kstr++)
        substr[ksubstr++] = str[kstr];
    substr[ksubstr] = '\0';
    return substr;
}

/*
    Finds the closing bracket in the null-terminated string \code{str}.

    Given a pair of matching brackets \code{open} and \code{close} and 
    an index $i$ such that \code{str[i] == open}, finds the index of the 
    corresponding closing bracket, or returns $-1$ if it doesn't exist.
 */
static size_t mpoly_str_find_close(const char * str, size_t i, char open, char close)
{
    size_t count = 1;

    while (str[++i] != '\0')
    {
        if (str[i] == close)
        {
            count--;
            if (count == 0)
                return i;
        }
        else if (str[i] == open)
            count++;
    }
    return -1;
}

/*
    Returns the monomial in $n$ variables given by the string 
    \code{str}, where the string is simply a space-separated 
    list of exponents.
 */
static mon_t mpoly_mon_set_str(const char * str, int n)
{
    mon_t rop;
    int i;
    size_t off;

    mon_init(rop);
    off = 0;
    for (i = 0; i < n; i++)
    {
        exp_t e = atoi(str + off);
        mon_set_exp(rop, i, e);
        off += 1 + (e >= 0) + (e >= 10) + (e >= 100);
    }

    return rop;
}

int mpoly_set_str(mpoly_t rop, const char * str, const ctx_t ctx)
{
    int i, j, n, num;
    const size_t len = strlen(str);

    /* Step 1.  Set the number of variables */

    n = atoi(str);

    mpoly_clear(rop, ctx);
    mpoly_init(rop, n, ctx);

    /* Step 2.  Count the number of terms */

    num = 1;
    for (i = 0; i < len; i++)
    {
        if (num > 0 && str[i] == '[')
            num = -num - 1;
        else if (num < 0 && str[i] == ']')
            num = -num;
    }
    num = num - 1;

    if (num < 0)
    {
        printf("ERROR (mpoly_set_str).  num = %d\n", num);
        abort();
    }
    if (num == 0)
        return 1;

    /* Step 3.  Process the terms one after the other */
    
    /* Read past the first integer, and skip following white space */
    j = 0;
    while (str[j] != ' ')
        j++;
    while (str[j] == ' ')
        j++;
    
    for (i = 0; i < num; i++)
    {
        char *s;
        int ins, jclose;
        mon_t m, m2;
        char *c;
        void *c2;

        /*
           First we set the coefficient and the monomial separately, moving 
           the index j just one past the closing bracket ']' of the monomial.
         */

        mon_init(m);
        c = malloc(ctx->size);
        ctx->init(c);

        while (str[j] != '(' && str[j] != '[')
            j++;

        if (str[j] == '(')
        {
            jclose = mpoly_str_find_close(str, j, '(', ')');
            s = mpoly_str_substr(str, j + 1, jclose);
            ctx->set_str(c, s);
            free(s);

            j = jclose + 1;
            while (str[j] != '[')
                j++;
        }
        else
        {
            ctx->one(c);
        }
        
        jclose = mpoly_str_find_close(str, j, '[', ']');
        s = mpoly_str_substr(str, j + 1, jclose);
        m = mpoly_mon_set_str(s, n);
        free(s);

        j = jclose + 1;

        /* Now we insert the new node */

        ins = RBTREE_INSERT(mpoly, &m2, &c2, rop->dict, m, c, &mon_cmp);

        if (ins)
        {
            printf("ERROR (mpoly_set_str).  Duplicate monomial.\n");
            abort();
        }
    }
    return 1;
}
