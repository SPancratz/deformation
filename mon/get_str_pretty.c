#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mon.h"

#define MON_VARS(n)\
    (n == 1 ? "X" : \
     n == 2 ? "XY" : \
     n == 3 ? "XYZ" : \
     n == 4 ? "WXYZ" : \
     n == 5 ? "VWXYZ" : \
     n == 6 ? "UVWXYZ" : \
     n == 7 ? "TUVWXYZ" : \
              "STUVWXYZ")

char * mon_get_str_pretty(mon_t x, int n, const char * vars)
{
    char * str;
    int i;
    size_t off;

    if (mon_is_one(x))
    {
        str = malloc(2);
        if (!str)
        {
            printf("ERROR (mon_get_str_pretty).  Memory allocation failed.\n");
            abort();
        }
        
        str[0] = '1';
        str[1] = '\0';
        return str;
    }

    str = malloc(n * 6);
    if (!str)
    {
        printf("ERROR (mon_get_str_pretty).  Memory allocation failed.\n");
        abort();
    }
    
    off = 0;
    for (i = 0; x; i++, x >>= MON_BITS_PER_EXP)
    {
        exp_t e = x & MON_BITMASK_BLOCK;

        if (e == 1)
        {
            str[off++] = vars ? vars[i] : MON_VARS(n)[i];
            str[off++] = ' ';
        }
        else if (e > 1)
        {
            str[off++] = vars ? vars[i] : MON_VARS(n)[i];
            str[off++] = '^';
            sprintf(str + off, "%lu", e);
            off += mon_exp_len(e);
            str[off++] = ' ';
        }
    }
    
    str[off-1] = '\0';
    return str;
}

