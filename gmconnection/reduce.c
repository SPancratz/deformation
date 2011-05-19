#include "gmconnection.h" 

/*
    Reduces the element $Q \Omega / P^k$ in de Rham cohomology.

    Assumes that $Q$ is a polynomial of degree \f$kd - (n+1)\f$ that is to 
    be reduced, where $k = 2, \dotsc, n, n+1$ is sufficiently large for 
    this to make sense.  In fact, we assume $\ell \leq k \leq u+1$.

    This method sets the elements array $R$ as $g_{\ell}, \dotsc, g_u$ 
    such that $g_i$ is a basis element in $B_i$.  The entries of $R$ 
    for $0, \dotsc, \ell-1$ are not touched.

    Note that the polynomial \c Q is modified.

    \param[out] R     Array of length \f$u+2\f$ for the reduction of the 
                      element in de Rham cohomology
    \param[in] Q      Numerator of the element in de Rham cohomology
    \param[in] k      Pole order of the element in de Rham cohomology
    \param[in] d      Degree of \f$P\f$
    \param[in] dP     Array of partial derivatives of \f$P\f$
    \param[in] AUX    Array of auxiliary matrices
    \param[in] AUX_b  Temporary work space
    \param[in] AUX_x  Temporary work space
    \param[in] AUX_y  Temporary work space
    \param[in] l      Lower bound \f$\ell\f$
    \param[in] u      Upper bound \f$u\f$
 */
void gmc_reduce(mpoly_t *R, 
                const mpoly_t Q0, long k, long d, mpoly_t *dP, 
                mat_csr_solve_t *s, 
                mon_t **rows, mon_t **cols, long **p, 
                long l, long u, 
                const ctx_t ctx)
{
    long i;
    long n;
    mpoly_t *A;
    mpoly_t dAdX;
    mpoly_t Q;

    /* Init */
    n = Q0->n;
    mpoly_init(Q, n, ctx);
    mpoly_set(Q, Q0, ctx);
    A = malloc(n * sizeof(mpoly_t));
    for (i = 0; i < n; i++)
        mpoly_init(A[i], n, ctx);
    mpoly_init(dAdX, n, ctx);
    
    /*
        Note that $\g_{n+1}$ is necessarily zero.
        Thus we first reduce to the case where k is at most n.
     */
    if (k == u + 1)
    {
        gmc_decompose_poly(A, Q, s[k], rows[k], cols[k], p[k], ctx);
        
        /* Set up the next polynomial to be decomposed */
        mpoly_zero(Q, ctx);
        for (i = 0; i < n; i++)
        {
            mpoly_derivative(dAdX, A[i], i, ctx);
            mpoly_add(Q, Q, dAdX, ctx);
        }
        k--;
        mpoly_scalar_div_si(Q, Q, k, ctx);
    }
    
    /* Ensure that all higher parts of the reduction array are zero */
    for (i = k + 1; i <= u; i++)
        mpoly_zero(R[i], ctx);
    
    while (!gmc_basis_contains(Q, d))
    {
        gmc_decompose_poly(A, Q, s[k], rows[k], cols[k], p[k], ctx);
        
        for (i = 0; i < n; i++)
            mpoly_submul(Q, A[i], dP[i], ctx);
        mpoly_swap(R[k], Q, ctx);
        
        /* Set up the next polynomial to be decomposed */
        mpoly_zero(Q, ctx);
        for (i = 0; i < n; i++)
        {
            mpoly_derivative(dAdX, A[i], i, ctx);
            mpoly_add(Q, Q, dAdX, ctx);
        }
        k--;
        mpoly_scalar_div_si(Q, Q, k, ctx);
    }
    
    /* Set the last element Q, which we know lies in the basis */
    if (!mpoly_is_zero(Q, ctx))
        mpoly_swap(R[k], Q, ctx);
    else if (k >= l)
        mpoly_zero(R[k], ctx);
    
    while (k > l)
        mpoly_zero(R[--k], ctx);
    
    /* Clear */
    for (i = 0; i < n; i++)
        mpoly_clear(A[i], ctx);
    free(A);
    mpoly_clear(dAdX, ctx);
    mpoly_clear(Q, ctx);
}

