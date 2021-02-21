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
    const uint cHuffmanMaxSupportedSyms = 8192;

    CRN_EXPORT void* create_generate_huffman_codes_tables();
    CRN_EXPORT void free_generate_huffman_codes_tables(void* p);

    CRN_EXPORT bool generate_huffman_codes(void* pContext, uint num_syms, const uint16* pFreq, uint8* pCodesizes, uint& max_code_size, uint& total_freq_ret);
} // namespace crnlib
