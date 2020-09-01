// File: corpus_test.h
// See Copyright Notice and license at the end of inc/crnlib.h
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
