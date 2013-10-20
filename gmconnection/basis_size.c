/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

long gmc_basis_size(long n, long d)
{
    long a = d - 1, ans, i;

    ans = a;
    for (i = 1; i < n; i++)
        ans *= a;
    if (n % 2)
        ans++;
    else
        ans--;
    ans /= d;
    ans *= a;
    return ans;
}

