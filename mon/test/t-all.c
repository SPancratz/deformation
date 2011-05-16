#include <sys/types.h>
#include <time.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "mon.h"

int main(int argc, char *argv[])
{
    mon_t m, m1, m2, *list;
    long i, len;
    char *out, *out1, *out2;

    printf("t-all\n");
    printf("=====\n");
    fflush(stdout);

    /* INIT, FROM_STRING, CLEAR */
    
    mon_init(m);
    mon_set_str(m, "3  1 1 1");
    
    out = mon_get_str(m, 3);
    printf("m = {%s}\n", out);
    free(out);
    
    mon_clear(m);
    
    /* GENERATE_BY_DEGREE */

    printf("Generating all monomials in X, Y, Z of degree 6\n");
    list = mon_generate_by_degree(&len, 3, 6);
    for (i = 0; i < len; i++)
    {
        out = mon_get_str_pretty(list[i], 3, NULL);
        printf("    %s\n", out);
        free(out);
    }
    for (i = 0; i < len; i++)
        mon_clear(list[i]);
    free(list);
    
    printf("Generating all monomials in X, Y, Z of degree 6\n");
    list = mon_generate_by_degree_invlex(&len, 3, 6);
    for (i = 0; i < len; i++)
    {
        out = mon_get_str_pretty(list[i], 3, NULL);
        printf("    %lu - %s\n", list[i], out);
        free(out);
    }
    for (i = 0; i < len; i++)
        mon_clear(list[i]);
    free(list);
    
    /* IS_ONE, DEGREE */

    mon_init(m);
    mon_set_str(m, "3  2 0 1");
    out = mon_get_str_pretty(m, 3, "XYZ");
    printf("Is %s equal to one?  %d\n", out, mon_is_one(m));
    printf("Degree: %d\n", mon_degree(m));
    free(out);
    mon_clear(m);
    
    mon_init(m);
    mon_set_str(m, "3  0 0 0");
    out = mon_get_str_pretty(m, 3, "XYZ");
    printf("Is %s equal to one?  %d\n", out, mon_is_one(m));
    printf("Degree: %d\n", mon_degree(m));
    free(out);
    mon_clear(m);
    
    /* MUL, DIVIDES, INVLEX */
    
    mon_init(m);
    mon_init(m1);
    mon_init(m2);
    mon_set_str(m1, "3  2 0 1");
    mon_set_str(m2, "3  0 1 0");
    mon_mul(m, m1, m2);
    out = mon_get_str_pretty(m, 3, "XYZ");
    printf("Product is [2 0 1] [0 1 0] = %s\n", out);
    free(out);
    printf("Divides?  %d\n", mon_divides(m1, m2));
    printf("INVLEX: %d\n", mon_cmp_invlex(m1, m2));
    mon_clear(m1);
    mon_clear(m2);
    mon_clear(m);
    
    /* DIVIDES, DIV, INVLEX */
    
    mon_init(m);
    mon_init(m1);
    mon_init(m2);
    mon_set_str(m1, "3  1 0 1");
    mon_set_str(m2, "3  1 1 3");
    mon_mul(m, m1, m2);
    out = mon_get_str_pretty(m, 3, "XYZ");
    printf("Product is [1 0 1] [1 1 3] = %s\n", out);
    free(out);
    printf("Divides?  %d\n", mon_divides(m1, m2));
    mon_div(m, m2, m1);
    out = mon_get_str_pretty(m, 3, "XYZ");
    printf("Quotient: %s\n", out);
    free(out);
    printf("INVLEX: %d\n", mon_cmp_invlex(m1, m2));
    mon_clear(m1);
    mon_clear(m2);
    mon_clear(m);
    
    /* INVLEX */
    
    mon_set_str(m1, "3  6 0 0");
    mon_set_str(m2, "3  0 0 6");
    out1 = mon_get_str_pretty(m1, 3, NULL);
    out2 = mon_get_str_pretty(m2, 3, NULL);
    printf("CMP_INVLEX(%s, %s) = %d\n", out1, out2, mon_cmp_invlex(m1, m2));
    free(out1);
    free(out2);

    printf("... ");

    printf("PASS\n");
    fflush(stdout);
    return EXIT_SUCCESS;
}

