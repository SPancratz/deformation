#include "diagfrob.h"

void diagfrob_coefficient_mpq(mpq_t rop, long m, long p)
{
    if (m == 0)
    {
        mpq_set_si(rop, 1, 1);
    }
    if (m == 1)
    {
        mpq_set_si(rop, p == 2 ? -2 : 1, 1);
    }
    else
    {
        const long low = (m + (p - 1)) / p;  /* Ceiling of m/p */
        long ell;   /* Index variable for factorials */
        long n, k;  /* Indices of the double sum */
        mpq_t r;    /* Term that is added to the coefficient at each step */
        mpq_t s;    /* Factor that changes the above term */
        
        mpq_init(r);
        mpq_init(s);
        mpq_set_si(rop, 0, 1);
        
        mpz_ui_pow_ui(mpq_numref(r), p, m / (p - 1));
        mpz_fac_ui(mpq_denref(r), m);
        mpq_canonicalize(r);
        mpq_set(rop, r);

        for (n = m - (p - 1), k = 1; n >= low; n -= (p - 1), k++)
        {
            mpz_set_si(mpq_numref(s), 1);
            for (ell = 1; ell <= p; ell++)
                mpz_mul_si(mpq_numref(s), mpq_numref(s), n - k + ell);
            mpz_set_si(mpq_denref(s), k);
            mpz_mul_si(mpq_denref(s), mpq_denref(s), p);
            mpq_canonicalize(s);
            mpq_mul(r, r, s);
            mpq_add(rop, rop, r);
        }

        if ((m / (p - 1)) & 1L)
            mpq_neg(rop, rop);
        
        mpq_clear(r);
        mpq_clear(s);
    }
}
