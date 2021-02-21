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

#include "crn_sysdetection.h"

#if defined(CRN_CC_MSVC)
#pragma warning(disable : 4201) // nonstandard extension used : nameless struct/union
#pragma warning(disable : 4127) // conditional expression is constant
#pragma warning(disable : 4793) // function compiled as native
#pragma warning(disable : 4324) // structure was padded due to __declspec(align())
#endif

#if defined(CRN_CC_MSVC)
#define CRN_NEVER_INLINE __declspec(noinline)
#define CRN_FORCE_INLINE __forceinline
#elif defined(CRN_CC_GNU)
#define CRN_NEVER_INLINE __attribute__((noinline))
#define CRN_FORCE_INLINE inline __attribute__((always_inline))
#else
#define CRN_NEVER_INLINE
#define CRN_FORCE_INLINE inline
#endif

#if defined(WIN32) && !defined(CRNLIB_ANSI_CPLUSPLUS)
// MSVC or MinGW, x86 or x64, Win32 API's for threading and Win32 Interlocked API's or GCC built-ins for atomic ops.
#ifdef NDEBUG
// Ensure checked iterators are disabled. Note: Be sure anything else that links against this lib also #define's this stuff, or remove this crap!
#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
#endif
#ifndef _DLL
// If we're using the DLL form of the run-time libs, we're also going to be enabling exceptions because we'll be building CLR apps.
// Otherwise, we disable exceptions for a small speed boost.
#define _HAS_EXCEPTIONS 0
#endif
#define NOMINMAX

#define CRNLIB_USE_WIN32_API 1

#if defined(__MINGW32__) || defined(__MINGW64__)
#define CRNLIB_USE_GCC_ATOMIC_BUILTINS 1
#else
#define CRNLIB_USE_WIN32_ATOMIC_FUNCTIONS 1
#endif

#define CRNLIB_PLATFORM_PC 1

#if defined(_WIN64) || defined(__MINGW64__) || defined(_LP64) || defined(__LP64__)
#define CRNLIB_PLATFORM_PC_X64 1
#define CRNLIB_64BIT_POINTERS 1
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 1
#define CRNLIB_LITTLE_ENDIAN_CPU 1
#else
#define CRNLIB_PLATFORM_PC_X86 1
#define CRNLIB_64BIT_POINTERS 0
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 0
#define CRNLIB_LITTLE_ENDIAN_CPU 1
#endif

#define CRNLIB_USE_UNALIGNED_INT_LOADS 1
#define CRNLIB_RESTRICT __restrict

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
#define CRNLIB_USE_MSVC_INTRINSICS 1
#endif

#define CRNLIB_INT64_FORMAT_SPECIFIER "%I64i"
#define CRNLIB_UINT64_FORMAT_SPECIFIER "%I64u"

#define CRNLIB_STDCALL __stdcall
#define CRNLIB_MEMORY_IMPORT_BARRIER
#define CRNLIB_MEMORY_EXPORT_BARRIER
#elif defined(__GNUC__) && !defined(CRNLIB_ANSI_CPLUSPLUS)
// GCC x86 or x64, pthreads for threading and GCC built-ins for atomic ops.
#define CRNLIB_PLATFORM_PC 1

#if defined(_WIN64) || defined(__MINGW64__) || defined(_LP64) || defined(__LP64__)
#define CRNLIB_PLATFORM_PC_X64 1
#define CRNLIB_64BIT_POINTERS 1
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 1
#else
#define CRNLIB_PLATFORM_PC_X86 1
#define CRNLIB_64BIT_POINTERS 0
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 0
#endif

#define CRNLIB_USE_UNALIGNED_INT_LOADS 1

#define CRNLIB_LITTLE_ENDIAN_CPU 1

#define CRNLIB_USE_PTHREADS_API 1
#define CRNLIB_USE_GCC_ATOMIC_BUILTINS 1

#define CRNLIB_RESTRICT

#define CRNLIB_INT64_FORMAT_SPECIFIER "%lli"
#define CRNLIB_UINT64_FORMAT_SPECIFIER "%llu"

#define CRNLIB_STDCALL
#define CRNLIB_MEMORY_IMPORT_BARRIER
#define CRNLIB_MEMORY_EXPORT_BARRIER
#else
// Vanilla ANSI-C/C++
// No threading support, unaligned loads are NOT okay.
#if defined(_WIN64) || defined(__MINGW64__) || defined(_LP64) || defined(__LP64__)
#define CRNLIB_64BIT_POINTERS 1
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 1
#else
#define CRNLIB_64BIT_POINTERS 0
#define CRNLIB_CPU_HAS_64BIT_REGISTERS 0
#endif

#define CRNLIB_USE_UNALIGNED_INT_LOADS 0

#if __BIG_ENDIAN__
#define CRNLIB_BIG_ENDIAN_CPU 1
#else
#define CRNLIB_LITTLE_ENDIAN_CPU 1
#endif

#define CRNLIB_USE_GCC_ATOMIC_BUILTINS 0
#define CRNLIB_USE_WIN32_ATOMIC_FUNCTIONS 0

#define CRNLIB_RESTRICT

#define CRNLIB_INT64_FORMAT_SPECIFIER "%I64i"
#define CRNLIB_UINT64_FORMAT_SPECIFIER "%I64u"

#define CRNLIB_STDCALL
#define CRNLIB_MEMORY_IMPORT_BARRIER
#define CRNLIB_MEMORY_EXPORT_BARRIER
#endif

#define CRNLIB_SLOW_STRING_LEN_CHECKS 1

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <algorithm>
#include <locale>
#include <memory.h>
#include <limits.h>
#include <algorithm>
#include <errno.h>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#define CRNLIB_FALSE (0)
#define CRNLIB_TRUE (1)
#define CRNLIB_MAX_PATH (260)

#ifdef _DEBUG
#define CRNLIB_BUILD_DEBUG
#else
#define CRNLIB_BUILD_RELEASE

#ifndef NDEBUG
#define NDEBUG
#endif

#ifdef DEBUG
#error DEBUG cannot be defined in CRNLIB_BUILD_RELEASE
#endif
#endif

#include "crn_types.h"
#include "crn_assert.h"
#include "crn_platform.h"
#include "crn_helpers.h"
#include "crn_traits.h"
#include "crn_mem.h"
#include "crn_math.h"
#include "crn_utils.h"
#include "crn_hash.h"
#include "crn_vector.h"
#include "crn_timer.h"
#include "crn_dynamic_string.h"
