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

namespace crnlib
{
    typedef unsigned long long timer_ticks;

    class CRN_EXPORT timer
    {
    public:
        timer();
        timer(timer_ticks start_ticks);

        void start();
        void start(timer_ticks start_ticks);

        void stop();

        double get_elapsed_secs() const;
        inline double get_elapsed_ms() const
        {
            return get_elapsed_secs() * 1000.0f;
        }
        timer_ticks get_elapsed_us() const;

        static void init();
        static inline timer_ticks get_ticks_per_sec()
        {
            return g_freq;
        }
        static timer_ticks get_init_ticks();
        static timer_ticks get_ticks();
        static double ticks_to_secs(timer_ticks ticks);
        static inline double ticks_to_ms(timer_ticks ticks)
        {
            return ticks_to_secs(ticks) * 1000.0f;
        }
        static inline double get_secs()
        {
            return ticks_to_secs(get_ticks());
        }
        static inline double get_ms()
        {
            return ticks_to_ms(get_ticks());
        }

    private:
        static timer_ticks g_init_ticks;
        static timer_ticks g_freq;
        static double g_inv_freq;

        timer_ticks m_start_time;
        timer_ticks m_stop_time;

        bool m_started : 1;
        bool m_stopped : 1;
    };

    // Prints object's lifetime to stdout
    class timed_scope
    {
    public:
        inline timed_scope(const char* pName = "timed_scope"):
            m_pName(pName)
        {
            m_tm.start();
        }

        inline double get_elapsed_secs() const
        {
            return m_tm.get_elapsed_secs();
        }
        inline double get_elapsed_ms() const
        {
            return m_tm.get_elapsed_ms();
        }

        const timer& get_timer() const
        {
            return m_tm;
        }
        timer& get_timer()
        {
            return m_tm;
        }

        inline ~timed_scope()
        {
            double secs = m_tm.get_elapsed_secs();
            printf("%s: %f secs, %f ms\n", m_pName, secs, secs * 1000.0f);
        }
    private:
        const char* m_pName;
        timer m_tm;
    };

}  // namespace crnlib
