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

#include "crnlib.h"
#include "crn_export.h"

namespace crnlib
{
    class mipmapped_texture;

    class itexture_comp
    {
        CRNLIB_NO_COPY_OR_ASSIGNMENT_OP(itexture_comp);

    public:
        itexture_comp()
        {
        }
        virtual ~itexture_comp()
        {
        }

        virtual const char* get_ext() const = 0;

        virtual bool compress_init(const crn_comp_params& params) = 0;
        virtual bool compress_pass(const crn_comp_params& params, float* pEffective_bitrate) = 0;
        virtual void compress_deinit() = 0;

        virtual const crnlib::vector<uint8>& get_comp_data() const = 0;
        virtual crnlib::vector<uint8>& get_comp_data() = 0;
    };

    CRN_EXPORT bool create_compressed_texture(const crn_comp_params& params, crnlib::vector<uint8>& comp_data,
        uint32* pActual_quality_level, float* pActual_bitrate);
    CRN_EXPORT bool create_texture_mipmaps(mipmapped_texture& work_tex, const crn_comp_params& params,
        const crn_mipmap_params& mipmap_params, bool generate_mipmaps);
    CRN_EXPORT bool create_compressed_texture(const crn_comp_params& params, const crn_mipmap_params& mipmap_params,
        crnlib::vector<uint8>& comp_data, uint32* pActual_quality_level,
        float* pActual_bitrate);
} // namespace crnlib
