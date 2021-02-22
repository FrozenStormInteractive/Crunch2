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
    struct chunk_tile_desc
    {
        // These values are in pixels, and always a multiple of cBlockPixelWidth/cBlockPixelHeight.
        uint m_x_ofs;
        uint m_y_ofs;
        uint m_width;
        uint m_height;
        uint m_layout_index;
    };

    struct chunk_encoding_desc
    {
        uint m_num_tiles;
        chunk_tile_desc m_tiles[4];
    };

    const uint cChunkPixelWidth = 8;
    const uint cChunkPixelHeight = 8;
    const uint cChunkBlockWidth = 2;
    const uint cChunkBlockHeight = 2;

    const uint cChunkMaxTiles = 4;

    const uint cBlockPixelWidthShift = 2;
    const uint cBlockPixelHeightShift = 2;

    const uint cBlockPixelWidth = 4;
    const uint cBlockPixelHeight = 4;

    const uint cNumChunkEncodings = 8;
    CRN_EXPORT extern chunk_encoding_desc g_chunk_encodings[cNumChunkEncodings];

    const uint cNumChunkTileLayouts = 9;
    const uint cFirst4x4ChunkTileLayout = 5;
    CRN_EXPORT extern chunk_tile_desc g_chunk_tile_layouts[cNumChunkTileLayouts];
} // namespace crnlib
