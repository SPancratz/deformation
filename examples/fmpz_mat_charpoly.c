/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "diagfrob.h"

int main(void)
{
    fmpz_mat_t A;
    fmpz_poly_t f;
    long n = 3;

    fmpz_mat_init(A, n, n);
    fmpz_poly_init(f);
    fmpz_t d;

    fmpz_set_ui(fmpz_mat_entry(A, 0, 0), 1);
    fmpz_set_ui(fmpz_mat_entry(A, 1, 1), 1);
    fmpz_set_ui(fmpz_mat_entry(A, 2, 2), 1);

    fmpz_mat_charpoly_modular(f, A);

    printf("A:\n");
    fmpz_mat_print_pretty(A);
    printf("\n");
    printf("f = ");
    fmpz_poly_print_pretty(f, "T");
    printf("\n");

    fmpz_set_si(fmpz_mat_entry(A, 0, 0),  2);
    fmpz_set_si(fmpz_mat_entry(A, 0, 1), -5);
    fmpz_set_si(fmpz_mat_entry(A, 0, 2),  1);
    fmpz_set_si(fmpz_mat_entry(A, 1, 0),  3);
    fmpz_set_si(fmpz_mat_entry(A, 1, 1), -1);
    fmpz_set_si(fmpz_mat_entry(A, 1, 2), -1);
    fmpz_set_si(fmpz_mat_entry(A, 2, 0), -1);
    fmpz_set_si(fmpz_mat_entry(A, 2, 1), -3);
    fmpz_set_si(fmpz_mat_entry(A, 2, 2),  0);

    fmpz_mat_charpoly_modular(f, A);

    printf("A:\n");
    fmpz_mat_print_pretty(A);
    printf("\n");
    printf("f = ");
    fmpz_poly_print_pretty(f, "T");
    printf("\n");

    fmpz_mat_clear(A);
    fmpz_poly_clear(f);

    return EXIT_SUCCESS;
}

