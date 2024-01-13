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

CRN_EXPORT bool crnlib_is_debugger_present(void);
CRN_EXPORT void crnlib_debug_break(void);
CRN_EXPORT void crnlib_output_debug_string(const char* p);

// actually in crnlib_assert.cpp
CRN_EXPORT void crnlib_assert(const char* pExp, const char* pFile, unsigned line);
CRN_EXPORT void crnlib_fail(const char* pExp, const char* pFile, unsigned line);

#if CRNLIB_LITTLE_ENDIAN_CPU
const bool c_crnlib_little_endian_platform = true;
#else
const bool c_crnlib_little_endian_platform = false;
#endif

const bool c_crnlib_big_endian_platform = !c_crnlib_little_endian_platform;

#if defined(CRN_OS_WIN) 
#define crn_fopen(pDstFile, f, m) fopen_s(pDstFile, f, m)
#define crn_fseek _fseeki64
#define crn_ftell _ftelli64
#elif defined(CRN_OS_LINUX)
#define crn_fopen(pDstFile, f, m) *(pDstFile) = fopen64(f, m)
#define crn_fseek fseeko64
#define crn_ftell ftello64
#else
#pragma message("Using fopen, ftell, fseek for file I/O - this may not support large files.")
#define crn_fopen(pDstFile, f, m) *(pDstFile) = fopen(f, m)
#define crn_fseek(s, o, w) fseek(s, static_cast<long>(o), w)
#define crn_ftell ftell
#endif

#if CRNLIB_USE_WIN32_API
#define CRNLIB_BREAKPOINT DebugBreak();
#define CRNLIB_BUILTIN_EXPECT(c, v) c
#elif defined(__GNUC__)
  #if defined(__aarch64__)
    #define CRNLIB_BREAKPOINT __builtin_trap();
  #else
    #define CRNLIB_BREAKPOINT asm("int $3");
  #endif
#define CRNLIB_BUILTIN_EXPECT(c, v) __builtin_expect(c, v)
#else
#define CRNLIB_BREAKPOINT
#define CRNLIB_BUILTIN_EXPECT(c, v) c
#endif

#if defined(__GNUC__)
#define CRNLIB_ALIGNED(x) __attribute__((aligned(x)))
#define CRNLIB_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define CRNLIB_ALIGNED(x) __declspec(align(x))
#define CRNLIB_NOINLINE __declspec(noinline)
#else
#define CRNLIB_ALIGNED(x)
#define CRNLIB_NOINLINE
#endif

#define CRNLIB_GET_ALIGNMENT(v) ((!sizeof(v)) ? 1 : (__alignof(v) ? __alignof(v) : sizeof(uint32)))

#ifndef _MSC_VER
int sprintf_s(char* buffer, size_t sizeOfBuffer, const char* format, ...);
int vsprintf_s(char* buffer, size_t sizeOfBuffer, const char* format, va_list args);
char* strlwr(char* p);
char* strupr(char* p);
#define _stricmp strcasecmp
#define _strnicmp strncasecmp
#endif

inline bool crnlib_is_little_endian()
{
    return c_crnlib_little_endian_platform;
}
inline bool crnlib_is_big_endian()
{
    return c_crnlib_big_endian_platform;
}

inline bool crnlib_is_pc()
{
#ifdef CRNLIB_PLATFORM_PC
    return true;
#else
    return false;
#endif
}

inline bool crnlib_is_x86()
{
#ifdef CRNLIB_PLATFORM_PC_X86
    return true;
#else
    return false;
#endif
}

inline bool crnlib_is_x64()
{
#ifdef CRNLIB_PLATFORM_PC_X64
    return true;
#else
    return false;
#endif
}
