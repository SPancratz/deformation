/*
    Copyright (C) 2010, 2011, 2012, 2013 Sebastian Pancratz
 */

#ifndef FLINT_EX_H
#define FLINT_EX_H


void 
_fmpz_mod_poly_compose_smod(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p);


#endif

