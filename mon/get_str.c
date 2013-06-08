/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mon.h"

char * mon_get_str(mon_t x, int n)
{
    char * str;
    int i;
    
    str = malloc(3 + n * 4 + 1);
    if (!str)
    {
        printf("ERROR (mon_get_str).  Memory allocation failed.\n");
        abort();
    }
    
    sprintf(str, "%d ", n);
    for (i = 0; i < n; i++)
    {
        exp_t e;
        size_t off;

        e = x & MON_BITMASK_BLOCK;
        off = strlen(str);
        sprintf(str + off, " %lu", e);
        x >>= MON_BITS_PER_EXP;
    }
    
    return str;
}
