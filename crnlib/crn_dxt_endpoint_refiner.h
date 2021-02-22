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

#include "crn_dxt.h"

namespace crnlib
{
    // TODO: Experimental/Not fully implemented
    class dxt_endpoint_refiner
    {
    public:
        dxt_endpoint_refiner();

        struct params
        {
            params() :
                m_block_index(0),
                m_pPixels(nullptr),
                m_num_pixels(0),
                m_pSelectors(nullptr),
                m_alpha_comp_index(0),
                m_error_to_beat(cUINT64_MAX),
                m_dxt1_selectors(true),
                m_perceptual(true),
                m_highest_quality(true)
            {
            }

            uint m_block_index;

            const color_quad_u8* m_pPixels;
            uint m_num_pixels;

            const uint8* m_pSelectors;

            uint m_alpha_comp_index;

            uint64 m_error_to_beat;

            bool m_dxt1_selectors;
            bool m_perceptual;
            bool m_highest_quality;
        };

        struct results
        {
            uint16 m_low_color;
            uint16 m_high_color;
            uint64 m_error;
        };

        bool refine(const params& p, results& r);

    private:
        const params* m_pParams;
        results* m_pResults;

        void optimize_dxt1(vec3F low_color, vec3F high_color);
        void optimize_dxt5(vec3F low_color, vec3F high_color);
    };
} // namespace crnlib
