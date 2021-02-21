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

#include "crn_vec.h"
#include "crn_export.h"

namespace crnlib
{
    template <typename T>
    class CRN_EXPORT vec_interval
    {
    public:
        enum { N = T::num_elements };
        typedef typename T::scalar_type scalar_type;

        inline vec_interval(const T& v)
        {
            m_bounds[0] = v;
            m_bounds[1] = v;
        }

        inline vec_interval(const T& low, const T& high)
        {
            m_bounds[0] = low;
            m_bounds[1] = high;
        }

        inline void clear()
        {
            m_bounds[0].clear();
            m_bounds[1].clear();
        }

        inline const T& operator[](uint i) const
        {
            CRNLIB_ASSERT(i < 2);
            return m_bounds[i];
        }

        inline T& operator[](uint i)
        {
            CRNLIB_ASSERT(i < 2);
            return m_bounds[i];
        }

    private:
        T m_bounds[2];
    };

    typedef vec_interval<vec1F> vec_interval1F;
    typedef vec_interval<vec2F> vec_interval2F;
    typedef vec_interval<vec3F> vec_interval3F;
    typedef vec_interval<vec4F> vec_interval4F;

    typedef vec_interval2F aabb2F;
    typedef vec_interval3F aabb3F;

}  // namespace crnlib
