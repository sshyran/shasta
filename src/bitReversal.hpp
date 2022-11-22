#ifndef SHASTA_BIT_REVERSAL_HPP
#define SHASTA_BIT_REVERSAL_HPP

// See https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel

#include "cstdint.hpp"

namespace shasta {
    inline uint16_t bitReversal(uint16_t);
    inline uint32_t bitReversal(uint32_t);
    inline uint64_t bitReversal(uint64_t);
}



inline uint16_t shasta::bitReversal(uint16_t x)
{
    const uint16_t m1 = uint16_t(0x5555);
    const uint16_t m2 = uint16_t(0x3333);
    const uint16_t m4 = uint16_t(0x0F0F);

    x = ((x >> 1) & m1) | ((x & m1) << 1);
    x = ((x >> 2) & m2) | ((x & m2) << 2);
    x = ((x >> 4) & m4) | ((x & m4) << 4);
    x = (x >> 8) | (x << 8);
    return x;
}



inline uint32_t shasta::bitReversal(uint32_t x)
{
    x = ((x >> 1) & 0x55555555) | ((x & 0x55555555) << 1);
    x = ((x >> 2) & 0x33333333) | ((x & 0x33333333) << 2);
    x = ((x >> 4) & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    x = ((x >> 8) & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    x = ( x >> 16) | ( x << 16);
    return x;
}



inline uint64_t shasta::bitReversal(uint64_t x)
{
    x = ((x >> 1)  & 0x5555555555555555UL) | ((x & 0x5555555555555555UL) << 1 );
    x = ((x >> 2)  & 0x3333333333333333UL) | ((x & 0x3333333333333333UL) << 2 );
    x = ((x >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((x & 0x0F0F0F0F0F0F0F0FUL) << 4 );
    x = ((x >> 8)  & 0x00FF00FF00FF00FFUL) | ((x & 0x00FF00FF00FF00FFUL) << 8 );
    x = ((x >> 16) & 0x0000FFFF0000FFFFUL) | ((x & 0x0000FFFF0000FFFFUL) << 16);
    x = (x >> 32) | (x << 32);
    return x;
}

#endif
