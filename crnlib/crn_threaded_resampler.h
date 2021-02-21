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

#include "crn_resampler.h"
#include "crn_vec.h"
#include "crn_export.h"

namespace crnlib
{
    class task_pool;
    class CRN_EXPORT threaded_resampler
    {
        CRNLIB_NO_COPY_OR_ASSIGNMENT_OP(threaded_resampler);

    public:
        threaded_resampler(task_pool& tp);
        ~threaded_resampler();

        enum pixel_format
        {
            cPF_Y_F32,
            cPF_RGBX_F32,
            cPF_RGBA_F32,

            cPF_Total
        };

        struct params
        {
            params()
            {
                clear();
            }

            void clear()
            {
                utils::zero_object(*this);

                m_boundary_op = Resampler::BOUNDARY_CLAMP;
                m_sample_low = 0.0f;
                m_sample_high = 255.0f;
                m_Pfilter_name = CRNLIB_RESAMPLER_DEFAULT_FILTER;
                m_filter_x_scale = 1.0f;
                m_filter_y_scale = 1.0f;
            }

            pixel_format m_fmt;

            const void* m_pSrc_pixels;
            uint m_src_width;
            uint m_src_height;
            uint m_src_pitch;

            void* m_pDst_pixels;
            uint m_dst_width;
            uint m_dst_height;
            uint m_dst_pitch;

            Resampler::Boundary_Op m_boundary_op;

            float m_sample_low;
            float m_sample_high;

            const char* m_Pfilter_name;
            float m_filter_x_scale;
            float m_filter_y_scale;
        };

        bool resample(const params& p);

    private:
        task_pool* m_pTask_pool;

        const params* m_pParams;

        Resampler::Contrib_List* m_pX_contribs;
        Resampler::Contrib_List* m_pY_contribs;
        uint m_bytes_per_pixel;

        crnlib::vector<vec4F> m_tmp_img;

        void free_contrib_lists();

        void resample_x_task(uint64 data, void* pData_ptr);
        void resample_y_task(uint64 data, void* pData_ptr);
    };
} // namespace crnlib
