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
#include "crn_utils.h"

namespace crnlib
{
    namespace utils
    {
        void endian_switch_words(uint16* p, uint num)
        {
            uint16* p_end = p + num;
            while (p != p_end)
            {
                uint16 k = *p;
                *p++ = swap16(k);
            }
        }

        void endian_switch_dwords(uint32* p, uint num)
        {
            uint32* p_end = p + num;
            while (p != p_end)
            {
                uint32 k = *p;
                *p++ = swap32(k);
            }
        }

        void copy_words(uint16* pDst, const uint16* pSrc, uint num, bool endian_switch)
        {
            if (!endian_switch)
            {
                memcpy(pDst, pSrc, num << 1U);
            }
            else
            {
                uint16* pDst_end = pDst + num;
                while (pDst != pDst_end)
                {
                    *pDst++ = swap16(*pSrc++);
                }
            }
        }

        void copy_dwords(uint32* pDst, const uint32* pSrc, uint num, bool endian_switch)
        {
            if (!endian_switch)
            {
                memcpy(pDst, pSrc, num << 2U);
            }
            else
            {
                uint32* pDst_end = pDst + num;
                while (pDst != pDst_end)
                {
                    *pDst++ = swap32(*pSrc++);
                }
            }
        }

        uint compute_max_mips(uint width, uint height)
        {
            if ((width | height) == 0)
            {
                return 0;
            }

            uint num_mips = 1;

            while ((width > 1U) || (height > 1U))
            {
                width >>= 1U;
                height >>= 1U;
                num_mips++;
            }

            return num_mips;
        }
    } // namespace utils
} // namespace crnlib
