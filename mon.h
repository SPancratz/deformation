/******************************************************************************

    Copyright (C) 2009, 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef MON_H
#define MON_H

typedef unsigned long mon_t;
typedef unsigned long exp_t;

#define MON_BITS_PER_EXP  8
#define MON_BITMASK_BLOCK  (0xFFul)
#define MON_BITMASK_ALL  (~ (0ul))

/* Initialization ************************************************************/

#define mon_init(x)  ((x) = 0)
#define mon_clear(x)
#define mon_set(x, y)  ((x) = (y))
#define mon_swap(x, y)  \
    do { mon_t xxx = (x); (x) = (y); (y) = xxx; } while(0)
#define mon_one(x)  ((x) = 0)

/* Access ********************************************************************/

#define mon_get_exp(x, i)  \
    ((exp_t) (((x) >> ((i) * MON_BITS_PER_EXP)) & MON_BITMASK_BLOCK))

#define mon_set_exp(x, i, e)  \
    ((x) = (((x) & (MON_BITMASK_ALL - (MON_BITMASK_BLOCK << ((i) * MON_BITS_PER_EXP)))) | (((mon_t)(e)) << ((i) * MON_BITS_PER_EXP))))

#define mon_inc_exp(x, i, e)  \
    ((x) = ((x) + (((mon_t)(e)) << ((i) * MON_BITS_PER_EXP))))

#define mon_dec_exp(x, i, e)  \
    ((x) = ((x) - (((mon_t)(e)) << ((i) * MON_BITS_PER_EXP))))

/* Comparison ****************************************************************/

#define mon_cmp_invlex(x, y)  (((x) < (y)) ? -1 : ((x) > (y) ? 1 : 0))

#define mon_is_one(x)  ((x) == 0)

/* Multiplication and division ***********************************************/

#define mon_mul(x, y, z)  ((x) = (y) + (z))

#define mon_div(x, y, z)  ((x) = (y) - (z))

int mon_divides(mon_t x, mon_t y);

/* Monomial parameters *******************************************************/

int mon_degree(mon_t x);

/* Input and output **********************************************************/

static int mon_exp_len(exp_t e)
{
    if (e < 10)
        return 1;
    else if (e < 100)
        return 2;
    else 
        return 3;
}

char * mon_get_str(mon_t x, int n);

char * mon_get_str_pretty(mon_t x, int n, const char * vars);

mon_t mon_set_str(char * str);

/* Enumeration  **************************************************************/

mon_t * mon_generate_by_degree(long * len, int n, int d);

mon_t * mon_generate_by_degree_invlex(long * len, int n, int d);

/* Auxiliary functions *******************************************************/

static unsigned long mon_binom(unsigned long n, unsigned long k)
{
    unsigned long i, d, res;

    if (n < k)
        return 0;
    if (k > n / 2)
        k = n - k;
    if (k == 0)
        return 1;
    d = n - k;
    res = d + 1;
    for (i = 2; i <= k; i++)
    {
        res *= d + i;
        res /= i;
    }
    return res;
}

#endif

