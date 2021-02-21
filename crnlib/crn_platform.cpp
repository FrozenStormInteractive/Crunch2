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

#ifndef _MSC_VER
int sprintf_s(char* buffer, size_t sizeOfBuffer, const char* format, ...)
{
    if (!sizeOfBuffer)
    {
        return 0;
    }

    va_list args;
    va_start(args, format);
    int c = vsnprintf(buffer, sizeOfBuffer, format, args);
    va_end(args);

    buffer[sizeOfBuffer - 1] = '\0';

    if (c < 0)
    {
        return sizeOfBuffer - 1;
    }

    return CRNLIB_MIN(c, (int)sizeOfBuffer - 1);
}

int vsprintf_s(char* buffer, size_t sizeOfBuffer, const char* format, va_list args)
{
    if (!sizeOfBuffer)
    {
        return 0;
    }

    int c = vsnprintf(buffer, sizeOfBuffer, format, args);

    buffer[sizeOfBuffer - 1] = '\0';

    if (c < 0)
    {
        return sizeOfBuffer - 1;
    }

    return CRNLIB_MIN(c, (int)sizeOfBuffer - 1);
}

char* strlwr(char* p)
{
    char* q = p;
    while (*q)
    {
        char c = *q;
        *q++ = tolower(c);
    }
    return p;
}

char* strupr(char* p)
{
    char* q = p;
    while (*q)
    {
        char c = *q;
        *q++ = toupper(c);
    }
    return p;
}
#endif

void crnlib_debug_break(void)
{
    CRNLIB_BREAKPOINT
}

#if CRNLIB_USE_WIN32_API
#include "crn_winhdr.h"

bool crnlib_is_debugger_present(void)
{
    return IsDebuggerPresent() != 0;
}

void crnlib_output_debug_string(const char* p)
{
    OutputDebugStringA(p);
}
#else
bool crnlib_is_debugger_present(void)
{
    return false;
}

void crnlib_output_debug_string(const char* p)
{
    puts(p);
}
#endif  // CRNLIB_USE_WIN32_API
