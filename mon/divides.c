#include "mon.h"

int mon_divides(mon_t x, mon_t y)
{
    while (x != 0)
    {
        if ((x & MON_BITMASK_BLOCK) > (y & MON_BITMASK_BLOCK))
            return 0;
        x >>= MON_BITS_PER_EXP;
        y >>= MON_BITS_PER_EXP;
    }
    return 1;
}
