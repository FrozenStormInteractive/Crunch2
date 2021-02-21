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
#include "crn_colorized_console.h"
#ifdef CRNLIB_USE_WIN32_API
#include "crn_winhdr.h"
#endif

namespace crnlib
{
    void colorized_console::init()
    {
        console::init();
        console::add_console_output_func(console_output_func, nullptr);
    }

    void colorized_console::deinit()
    {
        console::remove_console_output_func(console_output_func);
        console::deinit();
    }

    void colorized_console::tick()
    {
    }

#ifdef CRNLIB_USE_WIN32_API
    bool colorized_console::console_output_func(eConsoleMessageType type, const char* pMsg, void*)
    {
        if (console::get_output_disabled())
        {
            return true;
        }

        HANDLE cons = GetStdHandle(STD_OUTPUT_HANDLE);

        DWORD attr = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;
        switch (type)
        {
        case cDebugConsoleMessage:
            attr = FOREGROUND_BLUE | FOREGROUND_INTENSITY;
            break;
        case cMessageConsoleMessage:
            attr = FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY;
            break;
        case cWarningConsoleMessage:
            attr = FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY;
            break;
        case cErrorConsoleMessage:
            attr = FOREGROUND_RED | FOREGROUND_INTENSITY;
            break;
        default:
            break;
        }

        if (INVALID_HANDLE_VALUE != cons)
        {
            SetConsoleTextAttribute(cons, (WORD)attr);
        }

        if ((console::get_prefixes()) && (console::get_at_beginning_of_line()))
        {
            switch (type)
            {
            case cDebugConsoleMessage:
                printf("Debug: %s", pMsg);
                break;
            case cWarningConsoleMessage:
                printf("Warning: %s", pMsg);
                break;
            case cErrorConsoleMessage:
                printf("Error: %s", pMsg);
                break;
            default:
                printf("%s", pMsg);
                break;
            }
        }
        else
        {
            printf("%s", pMsg);
        }

        if (console::get_crlf())
        {
            printf("\n");
        }

        if (INVALID_HANDLE_VALUE != cons)
        {
            SetConsoleTextAttribute(cons, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        }

        return true;
    }
#else
    bool colorized_console::console_output_func(eConsoleMessageType type, const char* pMsg, void*)
    {
        if (console::get_output_disabled())
        {
            return true;
        }

        if ((console::get_prefixes()) && (console::get_at_beginning_of_line()))
        {
            switch (type)
            {
            case cDebugConsoleMessage:
                printf("Debug: %s", pMsg);
                break;
            case cWarningConsoleMessage:
                printf("Warning: %s", pMsg);
                break;
            case cErrorConsoleMessage:
                printf("Error: %s", pMsg);
                break;
            default:
                printf("%s", pMsg);
                break;
            }
        }
        else
        {
            printf("%s", pMsg);
        }

        if (console::get_crlf())
        {
            printf("\n");
        }

        return true;
    }
#endif

} // namespace crnlib
