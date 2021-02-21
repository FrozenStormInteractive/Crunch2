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

#include "crn_core.h"

#if CRNLIB_USE_WIN32_API
#include "crn_winhdr.h"
#endif

static bool g_fail_exceptions;
static bool g_exit_on_failure = true;

void crnlib_enable_fail_exceptions(bool enabled)
{
    g_fail_exceptions = enabled;
}

void crnlib_assert(const char* pExp, const char* pFile, unsigned line)
{
    char buf[512];

    sprintf_s(buf, sizeof(buf), "%s(%u): Assertion failed: \"%s\"\n", pFile, line, pExp);

    crnlib_output_debug_string(buf);

    fputs(buf, stderr);

    if (crnlib_is_debugger_present())
    {
        crnlib_debug_break();
    }
}

void crnlib_fail(const char* pExp, const char* pFile, unsigned line)
{
    char buf[512];

    sprintf_s(buf, sizeof(buf), "%s(%u): Failure: \"%s\"\n", pFile, line, pExp);

    crnlib_output_debug_string(buf);

    fputs(buf, stderr);

    if (crnlib_is_debugger_present())
    {
        crnlib_debug_break();
    }

#if CRNLIB_USE_WIN32_API
    if (g_fail_exceptions)
    {
        RaiseException(CRNLIB_FAIL_EXCEPTION_CODE, 0, 0, nullptr);
    }
    else
#endif
    {
        if (g_exit_on_failure)
        {
            exit(EXIT_FAILURE);
        }
    }
}

void trace(const char* pFmt, va_list args)
{
    if (crnlib_is_debugger_present())
    {
        char buf[512];
        vsprintf_s(buf, sizeof(buf), pFmt, args);

        crnlib_output_debug_string(buf);
    }
};

void trace(const char* pFmt, ...)
{
    va_list args;
    va_start(args, pFmt);
    trace(pFmt, args);
    va_end(args);
};
