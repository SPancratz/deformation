#include "diagfrob.h"

void diagfrob_falling_fac_mpq(mpq_t rop, long u, long d, long r)
{
    if (r == 0)
    {
        mpq_set_si(rop, 1, 1);
    }
    else
    {
        long j;
        mpq_t x;

        mpq_set_si(rop, u, d);
        mpq_canonicalize(rop);
        mpq_init(x);
        mpq_set(x, rop);
        for (j = 1; j < r; j++)
        {
            mpz_add(mpq_numref(x), mpq_numref(x), mpq_denref(x));
            if (mpq_sgn(x) == 0)
            {
                mpq_set_si(rop, 0, 1);
                break;
            }
            mpq_mul(rop, rop, x);
        }
        mpq_clear(x);
    }
}
