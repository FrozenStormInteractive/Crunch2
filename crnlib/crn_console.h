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
#include "crn_dynamic_string.h"
#include "crn_export.h"

#ifdef WIN32
#include <tchar.h>
#include <conio.h>
#endif

#if defined(CRN_CC_GNU)
#include <termios.h>
#include <unistd.h>
#endif

namespace crnlib
{
    class dynamic_string;
    class data_stream;
    class mutex;

    enum eConsoleMessageType
    {
        cDebugConsoleMessage, // debugging messages
        cProgressConsoleMessage, // progress messages
        cInfoConsoleMessage, // ordinary messages
        cConsoleConsoleMessage, // user console output
        cMessageConsoleMessage, // high importance messages
        cWarningConsoleMessage, // warnings
        cErrorConsoleMessage, // errors

        cCMTTotal,
    };

    typedef bool (*console_output_func)(eConsoleMessageType type, const char* pMsg, void* pData);

    class console
    {
    public:
        CRN_EXPORT static void init();
        CRN_EXPORT static void deinit();

        static bool is_initialized()
        {
            return m_pMutex != nullptr;
        }

        CRN_EXPORT static void set_default_category(eConsoleMessageType category);
        CRN_EXPORT static eConsoleMessageType get_default_category();

        CRN_EXPORT static void add_console_output_func(console_output_func pFunc, void* pData);
        CRN_EXPORT static void remove_console_output_func(console_output_func pFunc);

        CRN_EXPORT static void printf(const char* p, ...);

        CRN_EXPORT static void vprintf(eConsoleMessageType type, const char* p, va_list args);
        CRN_EXPORT static void printf(eConsoleMessageType type, const char* p, ...);

        CRN_EXPORT static void cons(const char* p, ...);
        CRN_EXPORT static void debug(const char* p, ...);
        CRN_EXPORT static void progress(const char* p, ...);
        CRN_EXPORT static void info(const char* p, ...);
        CRN_EXPORT static void message(const char* p, ...);
        CRN_EXPORT static void warning(const char* p, ...);
        CRN_EXPORT static void error(const char* p, ...);

        // FIXME: All console state is currently global!
        CRN_EXPORT static void disable_prefixes();
        CRN_EXPORT static void enable_prefixes();
        static bool get_prefixes()
        {
            return m_prefixes;
        }
        static bool get_at_beginning_of_line()
        {
            return m_at_beginning_of_line;
        }

        CRN_EXPORT static void disable_crlf();
        CRN_EXPORT static void enable_crlf();
        static bool get_crlf()
        {
            return m_crlf;
        }

        static void disable_output()
        {
            m_output_disabled = true;
        }
        static void enable_output()
        {
            m_output_disabled = false;
        }
        static bool get_output_disabled()
        {
            return m_output_disabled;
        }

        static void set_log_stream(data_stream* pStream)
        {
            m_pLog_stream = pStream;
        }
        static data_stream* get_log_stream()
        {
            return m_pLog_stream;
        }

        static uint get_num_messages(eConsoleMessageType type)
        {
            return m_num_messages[type];
        }

    private:
        static eConsoleMessageType m_default_category;

        struct console_func
        {
            console_func(console_output_func func = nullptr, void* pData = nullptr) :
                m_func(func),
                m_pData(pData)
            {
            }

            console_output_func m_func;
            void* m_pData;
        };

        CRN_EXPORT static crnlib::vector<console_func> m_output_funcs;

        CRN_EXPORT static bool m_crlf, m_prefixes, m_output_disabled;

        CRN_EXPORT static data_stream* m_pLog_stream;

        CRN_EXPORT static mutex* m_pMutex;

        CRN_EXPORT static uint m_num_messages[cCMTTotal];

        CRN_EXPORT static bool m_at_beginning_of_line;
    };

#if defined(WIN32)
    inline int crn_getch()
    {
        return _getch();
    }
#elif defined(CRN_CC_GNU)
    inline int crn_getch()
    {
        struct termios oldt, newt;
        int ch;
        tcgetattr(STDIN_FILENO, &oldt);
        newt = oldt;
        newt.c_lflag &= ~(ICANON | ECHO);
        tcsetattr(STDIN_FILENO, TCSANOW, &newt);
        ch = getchar();
        tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
        return ch;
    }
#else
    inline int crn_getch()
    {
        printf("crn_getch: Unimplemented");
        return 0;
    }
#endif
} // namespace crnlib
