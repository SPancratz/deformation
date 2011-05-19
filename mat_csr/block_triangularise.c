#include "mat_csr.h"

long _mat_csr_block_triangularise(long *arp, long *b, long n, const long *j, 
                                  const long *p, const long *lenr, long *w)
{
    long dummy, i, i1, i2, ii, k;
    long icnt, isn, ist, ist1, iv, iw, lcnt, num, nnm1, stp;

    int _goto80;

    long *lowl, *numb, *prev;

    lowl = w;
    numb = w + n;
    prev = w + 2 * n;

    /* 
       Explanation of the work space

       -  ``ARP`` - Array of length ``N``;  ``ARP[i]`` is one less than the 
          number of unsearched edges leaving node ``i``
       -  ``IB`` - Array of length ``N``;  ``IB[i]`` is the position in the 
          ordering of the start of the ``i``th block.  ``IB[N-i]`` holds 
          the node number of the ``i``th node on the stack
       -  ``LOWL`` - Array of length ``N``;  ``LOWL[i]`` is the smallest 
          stack positin of any node to which a path from node ``i`` has been 
          found.  It is set to ``N`` when node ``i`` is removed from the 
          stack
       -  ``NUMB`` - Array of length ``N``;  ``NUMB[i]`` is the position of 
          node ``i`` in the stack if it is on it, is the permuted order of node
          ``i`` for those nodes whose final position has been found and is 
          ``-1`` otherwise
       -  ``PREV`` - Array of length ``N``;  ``PREV[i]`` is the node at the end 
          of the path when node ``i`` was placed on the stack
     */

    /* Initialisation of arrays */
    for (i = 0; i < n; i++)
    {
        numb[i] = -1;
        arp[i]  = lenr[i] - 2;
    }

    /* icnt: number of nodes whose positions in final ordering are fixed */
    icnt = 0;

    /* num: number of blocks that have been found */
    num  = 0;
    nnm1 = n + n - 1;

    /* Added variables to control flow */
    _goto80  = 0;

    for (isn = 0; isn < n; isn++)
    {
        /* Look for a starting node */
        if (numb[isn] != -1)
            continue;

        iv = isn;

        /* ist: number of nodes on the stack; it is the stack pointer */
        ist = 0;

        /* Put node "iv" at beginning of the stack */
        lowl[iv]  = 0;
        numb[iv]  = 0;
        b[n - 1] = iv;

        /* LABEL 80 */
        /* The body of this loop puts a new node on the stack or backtracks */
        for (dummy = 0; dummy < nnm1; dummy++)
        {
            _goto80 = 0;
            
            i1 = arp[iv] + 1;
            
            /* Have all edges leaving node "iv" been searched? */
            if (i1 >= 0)
            {
                i2 = p[iv] + lenr[iv];
                i1 = i2 - i1;

                /*
                    Look at edges leaving node "IV" until one enters 
                    a new node or all edges are exhausted
                 */
                for (ii = i1-1; ii < i2; ii++)
                {
                    iw = j[ii];
                    
                    /* Has node IW been on stack already? */
                    if (numb[iw] == -1)
                    {
                        /* Put new node on the stack */
                        arp[iv]  = i2 - ii - 3;
                        prev[iw] = iv;
                        iv  = iw;
                        ist = ist + 1;
                        lowl[iv] = ist;
                        numb[iv] = ist;
                        k = n - ist - 1;
                        b[k] = iv;

                        _goto80 = 1;
                        break;
                    }
                       
                    /* Update value of lowl[iv] if necessary */
                    if (lowl[iw] < lowl[iv])
                        lowl[iv] = lowl[iw];
                }
                
                if (_goto80)
                    continue;
                
                /* There are no more edges leaving node "iv" */
                arp[iv] = -2;
            }

            /* Is node "iv" the root of a block? */
            if (lowl[iv] >= numb[iv])
            {
                /* Order nodes in a block */
                num  = num + 1;
                ist1 = n - ist;
                lcnt = icnt + 1;
                
                /*
                    Peel block off the top of the stack starting at the top 
                    and working down to the root of the block
                 */
                for (stp = ist1 - 1; stp < n; stp++)
                {
                    iw = b[stp];
                    lowl[iw] = n;
                    numb[iw] = icnt;
                    icnt = icnt + 1;
                    if (iw == iv)
                        break;
                }

                ist = n - stp - 2;
                b[num - 1] = lcnt - 1;
                
                /* Are there any nodes left on the stack? */
                if (ist == -1)
                {
                    /* Have all nodes been ordered? */
                    if (icnt < n)
                        break;
                    else
                        goto s_exit;
                }
            }
            
            /* Backtrack to previous node on path */
            iw = iv;
            iv = prev[iv];
            
            /* Update value of lowl[iv] if necessary */
            if (lowl[iw] < lowl[iv])
                lowl[iv] = lowl[iw];
        }
    }    

  s_exit:

    /* Put permutation in the required form */
    for (i = 0; i < n; i++)
        arp[numb[i]] = i;
    return num;
}

long mat_csr_block_triangularise(long *pi, long *b, const mat_csr_t A, 
                                                    const ctx_t ctx)
{
    long blocks, *w;

    w = malloc(3 * A->m * sizeof(long));

    if (!w)
    {
        printf("ERROR (mat_csr_block_triangularise).\n\n");
        abort();
    }

    blocks = _mat_csr_block_triangularise(pi, b, A->m, A->j, 
                                          A->p, A->lenr, w);

    free(w);
    return blocks;
}

