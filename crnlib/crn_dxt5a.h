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
#include "crn_export.h"

namespace crnlib
{
    class CRN_EXPORT dxt5_endpoint_optimizer
    {
    public:
        dxt5_endpoint_optimizer();

        struct params
        {
            params() :
                m_block_index(0),
                m_pPixels(nullptr),
                m_num_pixels(0),
                m_comp_index(3),
                m_quality(cCRNDXTQualityUber),
                m_use_both_block_types(true)
            {
            }

            uint m_block_index;

            const color_quad_u8* m_pPixels;
            uint m_num_pixels;
            uint m_comp_index;

            crn_dxt_quality m_quality;

            bool m_use_both_block_types;
        };

        struct results
        {
            uint8* m_pSelectors;

            uint64 m_error;

            uint8 m_first_endpoint;
            uint8 m_second_endpoint;

            uint8 m_block_type; // 1 if 6-alpha, otherwise 8-alpha
            bool m_reordered;
        };

        bool compute(const params& p, results& r);

    private:
        const params* m_pParams;
        results* m_pResults;

        crnlib::vector<uint8> m_unique_values;
        crnlib::vector<uint> m_unique_value_weights;

        crnlib::vector<uint8> m_trial_selectors;
        crnlib::vector<uint8> m_best_selectors;
        int m_unique_value_map[256];

        sparse_bit_array m_flags;

        void evaluate_solution(uint low_endpoint, uint high_endpoint);
    };

} // namespace crnlib
