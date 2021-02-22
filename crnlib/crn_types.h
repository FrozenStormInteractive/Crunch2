/*
 * Copyright (c) 2010-2016 Richard Geldreich, Jr. and Binomial LLC
 * Copyright (c) 2020 FrozenStorm Interactive, Yoann Potinet
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * 
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation or credits
 *    is required.
 * 
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 
 * 3. This notice may not be removed or altered from any source distribution.
 */

#pragma once

#include "crn_core.h"

namespace crnlib
{
    typedef unsigned char uint8;
    typedef signed char int8;
    typedef unsigned short uint16;
    typedef signed short int16;
    typedef unsigned int uint32;
    typedef uint32 uint;
    typedef signed int int32;

#if defined(CRN_CC_MSVC)
    typedef unsigned __int64 uint64;
    typedef signed __int64 int64;
#else
    typedef unsigned long long uint64;
    typedef long long int64;
#endif

    const uint8 cUINT8_MIN = 0;
    const uint8 cUINT8_MAX = 0xFFU;
    const uint16 cUINT16_MIN = 0;
    const uint16 cUINT16_MAX = 0xFFFFU;
    const uint32 cUINT32_MIN = 0;
    const uint32 cUINT32_MAX = 0xFFFFFFFFU;
    const uint64 cUINT64_MIN = 0;
    const uint64 cUINT64_MAX = 0xFFFFFFFFFFFFFFFFULL;  //0xFFFFFFFFFFFFFFFFui64;

    const int8 cINT8_MIN = -128;
    const int8 cINT8_MAX = 127;
    const int16 cINT16_MIN = -32768;
    const int16 cINT16_MAX = 32767;
    const int32 cINT32_MIN = (-2147483647 - 1);
    const int32 cINT32_MAX = 2147483647;
    const int64 cINT64_MIN = (int64)0x8000000000000000ULL;  //(-9223372036854775807i64 - 1);
    const int64 cINT64_MAX = (int64)0x7FFFFFFFFFFFFFFFULL;  // 9223372036854775807i64;

#if CRNLIB_64BIT_POINTERS
    typedef uint64 uint_ptr;
    typedef uint64 uint32_ptr;
    typedef int64 signed_size_t;
    typedef uint64 ptr_bits_t;
#else
    typedef unsigned int uint_ptr;
    typedef unsigned int uint32_ptr;
    typedef signed int signed_size_t;
    typedef uint32 ptr_bits_t;
#endif

    enum eVarArg { cVarArg };
    enum eClear { cClear };
    enum eNoClamp { cNoClamp };
    enum { cInvalidIndex = -1 };

    const uint cIntBits = 32;

    struct empty_type {};

}  // namespace crnlib
