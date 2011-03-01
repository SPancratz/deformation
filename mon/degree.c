#include "mon.h"

int mon_degree(mon_t op)
{
	int deg = op & MON_BITMASK_BLOCK;
	while (op != 0)
		deg += (op >>= MON_BITS_PER_EXP) & MON_BITMASK_BLOCK;
	return deg;
}
