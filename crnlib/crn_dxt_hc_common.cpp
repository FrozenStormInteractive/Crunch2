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

#include "crn_core.h"
#include "crn_dxt_hc_common.h"

namespace crnlib
{
    chunk_encoding_desc g_chunk_encodings[cNumChunkEncodings] = {
        { 1, { { 0, 0, 8, 8, 0 } } },

        { 2, { { 0, 0, 8, 4, 1 }, { 0, 4, 8, 4, 2 } } },
        { 2, { { 0, 0, 4, 8, 3 }, { 4, 0, 4, 8, 4 } } },

        { 3, { { 0, 0, 8, 4, 1 }, { 0, 4, 4, 4, 7 }, { 4, 4, 4, 4, 8 } } },
        { 3, { { 0, 4, 8, 4, 2 }, { 0, 0, 4, 4, 5 }, { 4, 0, 4, 4, 6 } } },

        { 3, { { 0, 0, 4, 8, 3 }, { 4, 0, 4, 4, 6 }, { 4, 4, 4, 4, 8 } } },
        { 3, { { 4, 0, 4, 8, 4 }, { 0, 0, 4, 4, 5 }, { 0, 4, 4, 4, 7 } } },

        { 4, { { 0, 0, 4, 4, 5 }, { 4, 0, 4, 4, 6 }, { 0, 4, 4, 4, 7 }, { 4, 4, 4, 4, 8 } } }
    };

    chunk_tile_desc g_chunk_tile_layouts[cNumChunkTileLayouts] = {
        // 2x2
        { 0, 0, 8, 8, 0 },

        // 2x1
        { 0, 0, 8, 4, 1 },
        { 0, 4, 8, 4, 2 },

        // 1x2
        { 0, 0, 4, 8, 3 },
        { 4, 0, 4, 8, 4 },

        // 1x1
        { 0, 0, 4, 4, 5 },
        { 4, 0, 4, 4, 6 },
        { 0, 4, 4, 4, 7 },
        { 4, 4, 4, 4, 8 }
    };
} // namespace crnlib
