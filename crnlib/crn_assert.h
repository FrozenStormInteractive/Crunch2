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

#include "crn_export.h"

const unsigned int CRNLIB_FAIL_EXCEPTION_CODE = 256U;
CRN_EXPORT void crnlib_enable_fail_exceptions(bool enabled);

CRN_EXPORT void crnlib_assert(const char* pExp, const char* pFile, unsigned line);
CRN_EXPORT void crnlib_fail(const char* pExp, const char* pFile, unsigned line);

#ifdef NDEBUG
#define CRNLIB_ASSERT(x) ((void)0)
#undef CRNLIB_ASSERTS_ENABLED
#else
#define CRNLIB_ASSERT(_exp) (void)((!!(_exp)) || (crnlib_assert(#_exp, __FILE__, __LINE__), 0))
#define CRNLIB_ASSERTS_ENABLED
#endif

#define CRNLIB_VERIFY(_exp) (void)((!!(_exp)) || (crnlib_assert(#_exp, __FILE__, __LINE__), 0))

#define CRNLIB_FAIL(msg) \
    do \
    { \
        crnlib_fail(#msg, __FILE__, __LINE__); \
    } while (0)

#define CRNLIB_ASSERT_OPEN_RANGE(x, l, h) CRNLIB_ASSERT((x >= l) && (x < h))
#define CRNLIB_ASSERT_CLOSED_RANGE(x, l, h) CRNLIB_ASSERT((x >= l) && (x <= h))

CRN_EXPORT void trace(const char* pFmt, va_list args);
CRN_EXPORT void trace(const char* pFmt, ...);

// Borrowed from boost libraries.
template<bool x>
struct crnlib_assume_failure;
template<>
struct crnlib_assume_failure<true>
{
    enum
    {
        blah = 1
    };
};
template<int x>
struct crnlib_assume_try
{
};

#define CRNLIB_JOINER_FINAL(a, b) a##b
#define CRNLIB_JOINER(a, b) CRNLIB_JOINER_FINAL(a, b)
#define CRNLIB_JOIN(a, b) CRNLIB_JOINER(a, b)
#define CRNLIB_ASSUME(p) typedef crnlib_assume_try<sizeof(crnlib_assume_failure<(bool)(p)>)> CRNLIB_JOIN(crnlib_assume_typedef, __COUNTER__)

#ifdef NDEBUG
template<typename T>
inline T crnlib_assert_range(T i, T)
{
    return i;
}
template<typename T>
inline T crnlib_assert_range_incl(T i, T)
{
    return i;
}
#else
template<typename T>
inline T crnlib_assert_range(T i, T m)
{
    CRNLIB_ASSERT((i >= 0) && (i < m));
    return i;
}
template<typename T>
inline T crnlib_assert_range_incl(T i, T m)
{
    CRNLIB_ASSERT((i >= 0) && (i <= m));
    return i;
}
#endif
