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

#include "crn_core.h"

#ifdef WIN32
#define CRNLIB_PATH_SEPERATOR_CHAR '\\'
#else
#define CRNLIB_PATH_SEPERATOR_CHAR '/'
#endif

namespace crnlib
{
    CRN_EXPORT char* crn_strdup(const char* pStr);
    CRN_EXPORT int crn_stricmp(const char* p, const char* q);

    CRN_EXPORT char* strcpy_safe(char* pDst, uint dst_len, const char* pSrc);

    CRN_EXPORT bool int_to_string(int value, char* pDst, uint len);
    CRN_EXPORT bool uint_to_string(uint value, char* pDst, uint len);

    CRN_EXPORT bool string_to_int(const char*& pBuf, int& value);

    CRN_EXPORT bool string_to_uint(const char*& pBuf, uint& value);

    CRN_EXPORT bool string_to_int64(const char*& pBuf, int64& value);
    CRN_EXPORT bool string_to_uint64(const char*& pBuf, uint64& value);

    CRN_EXPORT bool string_to_bool(const char* p, bool& value);

    CRN_EXPORT bool string_to_float(const char*& p, float& value, uint round_digit = 512U);

    CRN_EXPORT bool string_to_double(const char*& p, double& value, uint round_digit = 512U);
    CRN_EXPORT bool string_to_double(const char*& p, const char* pEnd, double& value, uint round_digit = 512U);
} // namespace crnlib
