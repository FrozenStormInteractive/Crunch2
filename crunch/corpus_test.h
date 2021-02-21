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

#include "crn_command_line_params.h"
#include "crn_image_utils.h"

namespace crn
{
    class corpus_tester
    {
    public:
        corpus_tester();

        bool test(const char* pCmd_line);

    private:
        void print_comparative_metric_stats(const crnlib::command_line_params& params, const crnlib::vector<crnlib::image_utils::error_metrics>& stats1, const crnlib::vector<crnlib::image_utils::error_metrics>& stats2, crnlib::uint num_blocks_x, crnlib::uint num_blocks_y);
        void print_metric_stats(const crnlib::vector<crnlib::image_utils::error_metrics>& stats, crnlib::uint num_blocks_x, crnlib::uint num_blocks_y);

        crnlib::image_u8 m_bad_block_img;
        crnlib::uint m_next_bad_block_index;
        crnlib::uint m_total_bad_block_files;

        void flush_bad_blocks();
        void add_bad_block(crnlib::image_u8& block);
    };
} // namespace crnlib
