#include "gmconnection.h"

long gmc_basis_size(long n, long d)
{
    long a, ans, i;

    a = d - 1;
    ans = a;
    for (i = 2; i < n; i++)
        ans *= a;
    if (n % 2)
        ans--;
    else
        ans++;
    ans /= d;
    ans *= a;
    return ans;
}

