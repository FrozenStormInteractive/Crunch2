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
#include "crn_vector.h"

namespace crnlib
{
    bool elemental_vector::increase_capacity(uint min_new_capacity, bool grow_hint, uint element_size, object_mover pMover, bool nofail)
    {
        CRNLIB_ASSERT(m_size <= m_capacity);
#ifdef CRNLIB_64BIT_POINTERS
        CRNLIB_ASSERT(min_new_capacity < (0x400000000ULL / element_size));
#else
        CRNLIB_ASSERT(min_new_capacity < (0x7FFF0000U / element_size));
#endif

        if (m_capacity >= min_new_capacity)
        {
            return true;
        }

        ptr_bits_t new_capacity = min_new_capacity;
        if ((grow_hint) && (!math::is_power_of_2((uint64)new_capacity)))
        {
            new_capacity = math::next_pow2((uint64)new_capacity);
        }

        CRNLIB_ASSERT(new_capacity && (new_capacity > m_capacity));

        const size_t desired_size = element_size * new_capacity;
        size_t actual_size;
        if (!pMover)
        {
            void* new_p = crnlib_realloc(m_p, desired_size, &actual_size, true);
            if (!new_p)
            {
                if (nofail)
                {
                    return false;
                }

                char buf[256];
#ifdef _MSC_VER
                sprintf_s(buf, sizeof(buf), "vector: crnlib_realloc() failed allocating %u bytes", (uint)desired_size);
#else
                sprintf(buf, "vector: crnlib_realloc() failed allocating %u bytes", (uint)desired_size);
#endif
                CRNLIB_FAIL(buf);
            }
            m_p = new_p;
        }
        else
        {
            void* new_p = crnlib_malloc(desired_size, &actual_size);
            if (!new_p)
            {
                if (nofail)
                {
                    return false;
                }

                char buf[256];
#ifdef _MSC_VER
                sprintf_s(buf, sizeof(buf), "vector: crnlib_malloc() failed allocating %u bytes", (uint)desired_size);
#else
                sprintf(buf, "vector: crnlib_malloc() failed allocating %u bytes", (uint)desired_size);
#endif
                CRNLIB_FAIL(buf);
            }

            (*pMover)(new_p, m_p, m_size);

            if (m_p)
            {
                crnlib_free(m_p);
            }

            m_p = new_p;
        }

        if (actual_size > desired_size)
        {
            m_capacity = static_cast<uint>(actual_size / element_size);
        }
        else
        {
            m_capacity = static_cast<uint>(new_capacity);
        }

        return true;
    }
} // namespace crnlib
