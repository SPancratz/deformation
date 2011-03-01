#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mon.h"

mon_t mon_set_str(char * str)
{
    int i, n;   /* Number of variables */
    mon_t rop;  /* Output */
    exp_t e;    /* Temp */
    size_t off;

    n = atoi(str);
    if (!(1 <= n && n <= 8))
    {
        printf("ERROR (mon_set_str).  Bad number of variables.\n");
        abort();
    }
    
    mon_init(rop);
    off = 3;
    for (i = 0; i < n; i++)
    {
        exp_t e = atoi(str + off);

        mon_set_exp(rop, i, e);
        off += mon_exp_len(e) + 1;
    }
    
    return rop;
}
