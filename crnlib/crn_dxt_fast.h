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

#include "crn_color.h"
#include "crn_dxt.h"
#include "crn_export.h"

namespace crnlib
{
    namespace dxt_fast
    {
        CRN_EXPORT void compress_color_block(uint n, const color_quad_u8* block, uint& low16, uint& high16, uint8* pSelectors, bool refine = false);
        CRN_EXPORT void compress_color_block(dxt1_block* pDXT1_block, const color_quad_u8* pBlock, bool refine = false);

        CRN_EXPORT void compress_alpha_block(uint n, const color_quad_u8* block, uint& low8, uint& high8, uint8* pSelectors, uint comp_index);
        CRN_EXPORT void compress_alpha_block(dxt5_block* pDXT5_block, const color_quad_u8* pBlock, uint comp_index);

        CRN_EXPORT void find_representative_colors(uint n, const color_quad_u8* pBlock, color_quad_u8& lo, color_quad_u8& hi);
    } // namespace dxt_fast
} // namespace crnlib
