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

#include "crn_comp.h"
#include "crn_mipmapped_texture.h"
#include "crn_texture_comp.h"
#include "crn_export.h"

namespace crnlib
{
    class CRN_EXPORT dds_comp : public itexture_comp
    {
        CRNLIB_NO_COPY_OR_ASSIGNMENT_OP(dds_comp);

    public:
        dds_comp();
        virtual ~dds_comp();

        virtual const char* get_ext() const
        {
            return "DDS";
        }

        virtual bool compress_init(const crn_comp_params& params);
        virtual bool compress_pass(const crn_comp_params& params, float* pEffective_bitrate);
        virtual void compress_deinit();

        virtual const crnlib::vector<uint8>& get_comp_data() const
        {
            return m_comp_data;
        }
        virtual crnlib::vector<uint8>& get_comp_data()
        {
            return m_comp_data;
        }

    private:
        mipmapped_texture m_src_tex;
        mipmapped_texture m_packed_tex;

        crnlib::vector<uint8> m_comp_data;

        const crn_comp_params* m_pParams;

        pixel_format m_pixel_fmt;
        dxt_image::pack_params m_pack_params;

        task_pool m_task_pool;
        qdxt1_params m_q1_params;
        qdxt5_params m_q5_params;
        mipmapped_texture::qdxt_state* m_pQDXT_state;

        void clear();
        bool create_dds_tex(mipmapped_texture& dds_tex);
        bool convert_to_dxt(const crn_comp_params& params);
    };

} // namespace crnlib
