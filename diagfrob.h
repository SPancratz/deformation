/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef DIAGFROB_H
#define DIAGFROB_H

#include <stdlib.h>
#include <stdio.h>

#include "generics.h"
#include "flint.h"

#include "mat.h"

#include "mon.h"

#define DIAGFROB_MOD(x, m)  (((x) % (m) >= 0) ? (x) % (m) : (x) % (m) + (m))

