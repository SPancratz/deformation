#include "mat_csr.h"

long _mat_csr_zfdiagonal(long *pi, long n, const long *j, const long *p, 
                                   const long *lenr, long *w)
{
    long i1, i2, j1, j2, jord, k1, k2;
    long in1, in2;
    long numnz;

    long *pr, *arp, *cv, *out;

    /* Dummy initialisations, for compiler warnings */
    i2 = in2 = 0;

    /*
        Explanation of the work space

        -  ``pr[i]`` is the previous row to ``i`` in the depth first search.
        -  ``arp[i]`` is one less than the number of non-zero elements in 
           row ``i`` which have not been scanned when looking for a cheap 
           assignment
        -  ``cv[i]`` is the most recent row extension at which column ``i`` 
           was visited
        -  ``out[i]`` is one less than the number of non-zero elements in 
           row ``i`` which have not been scanned during one pass through the 
           main loop
     */

    pr  = w;
    arp = w + n;
    cv  = w + 2 * n;
    out = w + 3 * n;

    for (i1 = 0; i1 < n; i1++)
    {
        arp[i1] = lenr[i1] - 1;
        cv[i1]  = -1;
        pi[i1]  = -1;
    }

    numnz = 0;
    
    /*
        Main loop:  each pass round this loop either results 
        in a new assignment or gives a row with no assignment
     */
    for (jord = 0; jord < n; jord++)
    {
        j1 = jord;
        pr[j1] = -2;
        
        for (k1 = 0; k1 <= jord; k1++)
        {
            /* Look for a cheap assignment */
            in1 = arp[j1];
            if (in1 >= 0)
            {
                in2 = p[j1] + lenr[j1];
                in1 = in2 - in1;
                for (i2 = in1 - 1; i2 < in2; i2++)
                {
                    i1 = j[i2];
                    if (pi[i1] == -1)
                        goto s_110;
                }
                
                /* No cheap assignment in row */
                arp[j1] = -1;
            }

            /* Begin looking for assignment chain starting with row j1 */
            out[j1] = lenr[j1] - 1;

            /* Inner loop:  extends chain by one or backtracks */
            for (k2 = 0; k2 <= jord; k2++)
            {
                in1 = out[j1];

                if (in1 >= 0)
                {
                    in2 = p[j1] + lenr[j1];
                    in1 = in2 - in1;

                    /* Forward scan */
                    for (i2 = in1 - 1; i2 < in2; i2++)
                    {
                        i1 = j[i2];
                        if (cv[i1] == jord)
                            continue;

                        /* Column i1 not yet accessed during this pass */
                        j2 = j1;
                        j1 = pi[i1];
                        cv[i1]  = jord;
                        pr[j1]  = j2;
                        out[j2] = in2 - i2 - 2;

                        goto s_100;
                    }
                }
                
                /* Backtracking step */
                j1 = pr[j1];
                if (j1 == -2)
                    goto s_130;
            }

          s_100: ;
        }

      s_110:

        /* New assignment */
        pi[i1]  = j1;
        arp[j1] = in2 - i2 - 2;
        numnz = numnz + 1;
        for (k1 = 0; k1 <= jord; k1++)
        {
            j1 = pr[j1];
            if (j1 == -2)
                break;
            i2 = p[j1] + lenr[j1] - out[j1] - 2;
            i1 = j[i2];
            pi[i1] = j1;
        }

      s_130: ;
    }
    
    /* If the matrix is structurally singular, complete pi */
    if (numnz < n)
    {
        for (i1 = 0; i1 < n; i1++)
            arp[i1] = 0;
        k1 = 0;
        for (i1 = 0; i1 < n; i1++)
            if (pi[i1] == -1)
                out[k1++] = i1 + 1;
            else
                arp[pi[i1]] = i1 + 1;
        k1 = 0;
        for (i1 = 0; i1 < n; i1++)
            if (arp[i1] == 0)
                pi[out[k1++] - 1] = i1;
    }
    
    return numnz;
}

long mat_csr_zfdiagonal(long *pi, const mat_csr_t A)
{
    long numnz, *w;

    if (!(A->m == A->n))
    {
        printf("ERROR (zfdiagonal).\n\n");
        abort();
    }

    w = malloc(4 * A->m * sizeof(long));

    if (!w)
    {
        printf("ERROR (zfdiagonal).\n\n");
        abort();
    }

    numnz = _mat_csr_zfdiagonal(pi, A->m, A->j, A->p, A->lenr, w);

    free(w);

    return numnz;
}
