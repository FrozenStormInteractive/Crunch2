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
    const uint cInitAdler32 = 1U;
    CRN_EXPORT uint adler32(const void* pBuf, size_t buflen, uint adler32 = cInitAdler32);

    // crc16() intended for small buffers - doesn't use an acceleration table.
    const uint cInitCRC16 = 0;
    CRN_EXPORT uint16 crc16(const void* pBuf, size_t len, uint16 crc = cInitCRC16);

} // namespace crnlib
