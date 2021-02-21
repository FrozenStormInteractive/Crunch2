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

#include <LzmaLib.h>

#include "crn_core.h"
#include "crn_lzma_codec.h"
#include "crn_strutils.h"
#include "crn_checksum.h"
#include "crn_threading.h"

namespace crnlib
{
    lzma_codec::lzma_codec() :
        m_pCompress(LzmaCompress),
        m_pUncompress(LzmaUncompress)
    {
        CRNLIB_ASSUME(cLZMAPropsSize == LZMA_PROPS_SIZE);
    }

    lzma_codec::~lzma_codec()
    {
    }

    bool lzma_codec::pack(const void* p, uint n, crnlib::vector<uint8>& buf)
    {
        if (n > 1024U * 1024U * 1024U)
        {
            return false;
        }

        uint max_comp_size = n + math::maximum<uint>(128, n >> 8);
        buf.resize(sizeof(header) + max_comp_size);

        header* pHDR = reinterpret_cast<header*>(&buf[0]);
        uint8* pComp_data = &buf[sizeof(header)];

        utils::zero_object(*pHDR);

        pHDR->m_uncomp_size = n;
        pHDR->m_adler32 = adler32(p, n);

        if (n)
        {
            size_t destLen = 0;
            size_t outPropsSize = 0;
            int status = SZ_ERROR_INPUT_EOF;

            for (uint trial = 0; trial < 3; trial++)
            {
                destLen = max_comp_size;
                outPropsSize = cLZMAPropsSize;

                status = (*m_pCompress)(pComp_data, &destLen, reinterpret_cast<const unsigned char*>(p), n,
                    pHDR->m_lzma_props, &outPropsSize,
                    -1, /* 0 <= level <= 9, default = 5 */
                    0, /* default = (1 << 24) */
                    -1, /* 0 <= lc <= 8, default = 3  */
                    -1, /* 0 <= lp <= 4, default = 0  */
                    -1, /* 0 <= pb <= 4, default = 2  */
                    -1, /* 5 <= fb <= 273, default = 32 */
#ifdef WIN32
                    (g_number_of_processors > 1) ? 2 : 1
#else
                    1
#endif
                );

                if (status != SZ_ERROR_OUTPUT_EOF)
                {
                    break;
                }

                max_comp_size += ((n + 1) / 2);
                buf.resize(sizeof(header) + max_comp_size);
                pHDR = reinterpret_cast<header*>(&buf[0]);
                pComp_data = &buf[sizeof(header)];
            }

            if (status != SZ_OK)
            {
                buf.clear();
                return false;
            }

            pHDR->m_comp_size = static_cast<uint>(destLen);

            buf.resize(CRNLIB_SIZEOF_U32(header) + static_cast<uint32>(destLen));
        }

        pHDR->m_sig = header::cSig;
        pHDR->m_checksum = static_cast<uint8>(adler32((uint8*)pHDR + header::cChecksumSkipBytes, sizeof(header) - header::cChecksumSkipBytes));

        return true;
    }

    bool lzma_codec::unpack(const void* p, uint n, crnlib::vector<uint8>& buf)
    {
        buf.resize(0);

        if (n < sizeof(header))
        {
            return false;
        }

        const header& hdr = *static_cast<const header*>(p);
        if (hdr.m_sig != header::cSig)
        {
            return false;
        }

        if (static_cast<uint8>(adler32((const uint8*)&hdr + header::cChecksumSkipBytes, sizeof(hdr) - header::cChecksumSkipBytes)) != hdr.m_checksum)
        {
            return false;
        }

        if (!hdr.m_uncomp_size)
        {
            return true;
        }

        if (!hdr.m_comp_size)
        {
            return false;
        }

        if (hdr.m_uncomp_size > 1024U * 1024U * 1024U)
        {
            return false;
        }

        if (!buf.try_resize(hdr.m_uncomp_size))
        {
            return false;
        }

        const uint8* pComp_data = static_cast<const uint8*>(p) + sizeof(header);
        size_t srcLen = n - sizeof(header);
        if (srcLen < hdr.m_comp_size)
        {
            return false;
        }

        size_t destLen = hdr.m_uncomp_size;

        int status = (*m_pUncompress)(&buf[0], &destLen, pComp_data, &srcLen, hdr.m_lzma_props, cLZMAPropsSize);

        if ((status != SZ_OK) || (destLen != hdr.m_uncomp_size))
        {
            buf.clear();
            return false;
        }

        if (adler32(&buf[0], buf.size()) != hdr.m_adler32)
        {
            buf.clear();
            return false;
        }

        return true;
    }

} // namespace crnlib
