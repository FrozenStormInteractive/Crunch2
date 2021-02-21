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
    class CRN_EXPORT sparse_bit_array
    {
    public:
        sparse_bit_array();
        sparse_bit_array(uint size);
        sparse_bit_array(sparse_bit_array& other);
        ~sparse_bit_array();

        sparse_bit_array& operator=(sparse_bit_array& other);

        void clear();

        inline uint get_size()
        {
            return (m_num_groups << cBitsPerGroupShift);
        }

        void resize(uint size);

        sparse_bit_array& operator&=(const sparse_bit_array& other);
        sparse_bit_array& operator|=(const sparse_bit_array& other);
        sparse_bit_array& and_not(const sparse_bit_array& other);

        void swap(sparse_bit_array& other);

        void optimize();

        void set_bit_range(uint index, uint num);
        void clear_bit_range(uint index, uint num);

        void clear_all_bits();

        inline void set_bit(uint index)
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                pGroup = alloc_group(true);
                m_ppGroups[group_index] = pGroup;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            pGroup[bit_ofs >> 5] |= (1U << (bit_ofs & 31));
        }

        inline void clear_bit(uint index)
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                pGroup = alloc_group(true);
                m_ppGroups[group_index] = pGroup;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            pGroup[bit_ofs >> 5] &= (~(1U << (bit_ofs & 31)));
        }

        inline void set(uint index, bool value)
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                pGroup = alloc_group(true);
                m_ppGroups[group_index] = pGroup;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            uint bit = (1U << (bit_ofs & 31));

            uint c = pGroup[bit_ofs >> 5];
            uint mask = (uint)(-(int)value);

            pGroup[bit_ofs >> 5] = (c & ~bit) | (mask & bit);
        }

        inline bool get_bit(uint index) const
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                return 0;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            uint bit = (1U << (bit_ofs & 31));

            return (pGroup[bit_ofs >> 5] & bit) != 0;
        }

        inline uint32 get_uint32(uint index) const
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                return 0;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            return pGroup[bit_ofs >> 5];
        }

        inline void set_uint32(uint index, uint32 value) const
        {
            uint group_index = index >> cBitsPerGroupShift;
            CRNLIB_ASSERT(group_index < m_num_groups);

            uint32* pGroup = m_ppGroups[group_index];
            if (!pGroup)
            {
                pGroup = alloc_group(true);
                m_ppGroups[group_index] = pGroup;
            }

            uint bit_ofs = index & (cBitsPerGroup - 1);

            pGroup[bit_ofs >> 5] = value;
        }

        int find_first_set_bit(uint index, uint num) const;

        enum
        {
            cDWORDsPerGroupShift = 4U,
            cDWORDsPerGroup = 1U << cDWORDsPerGroupShift,

            cBitsPerGroupShift = cDWORDsPerGroupShift + 5,
            cBitsPerGroup = 1U << cBitsPerGroupShift,
            cBitsPerGroupMask = cBitsPerGroup - 1U,

            cBytesPerGroup = cDWORDsPerGroup * sizeof(uint32)
        };

        uint get_num_groups() const
        {
            return m_num_groups;
        }
        uint32** get_groups()
        {
            return m_ppGroups;
        }

    private:
        uint m_num_groups;
        uint32** m_ppGroups;

        static inline uint32* alloc_group(bool clear)
        {
            uint32* p = (uint32*)crnlib_malloc(cBytesPerGroup);
            CRNLIB_VERIFY(p);
            if (clear)
            {
                memset(p, 0, cBytesPerGroup);
            }
            return p;
        }

        static inline void free_group(void* p)
        {
            if (p)
            {
                crnlib_free(p);
            }
        }
    };

} // namespace crnlib
