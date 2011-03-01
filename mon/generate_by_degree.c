#include <stdlib.h>
#include <stdio.h>

#include "mon.h"

mon_t * mon_generate_by_degree(long * len, int n, int d)
{
    mon_t * st;  /* Stack (from the bottom), list (from the top) */
    int * sums;
    int * inds;
    int i;       /* Pointer to the top of the stack */
    int j;       /* Pointer to the bottom of the list */
    int N;       /* $\binom{n-1+d}{d}$ */
    int sum, ind;
    mon_t m;
    exp_t exp;
    
    /* Special case d == 0 */
    if (d == 0)
    {
        st = malloc(sizeof(mon_t));
        if (!st)
        {
            printf("ERROR (mon_generate_by_degree).  Memory allocation failed.\n");
            abort();
        }
        
        mon_init(st[0]);
        *len = 1;
        return st;
    }
    
    N = mon_binom(n - 1 + d, d);
    
    st   = malloc(N * sizeof(mon_t));
    sums = malloc(N * sizeof(int));
    inds = malloc(N * sizeof(int));
    if (!st || !sums || !inds)
    {
        printf("ERROR (mon_generate_by_degree).  Memory allocation failed.\n");
        abort();
    }
    
    for (i = 0; i < N; i++)
        mon_init(st[i]);
    
    mon_init(m);
    
    i = 0;
    sums[0] = 0;
    inds[0] = 0;
    j = N;
    
    while (i >= 0)
    {
        mon_set(m, st[i]);
        ind = inds[i];
        sum = sums[i];
        i--;
        
        if (ind < n - 1)
        {
            for (exp = 0; exp <= d - sum; exp++)
            {
                i++;
                mon_set(st[i], m);
                mon_inc_exp(st[i], ind, exp);
                inds[i] = ind + 1;
                sums[i] = sum + exp;
            }
        }
        else
        {
            j--;
            mon_set(st[j], m);
            mon_inc_exp(st[j], ind, d-sum);
        }
    }
    
    mon_clear(m);
    free(sums);
    free(inds);
    
    if (N > 1)
    {
        for (i = 0, j = N - 1; i < N / 2; i++, j--)
            mon_swap(st[i], st[j]);
    }
    
    *len = N;
    return st;
}
