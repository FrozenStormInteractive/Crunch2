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

namespace crnlib
{
    template<typename vector_type>
    class ray
    {
    public:
        typedef vector_type vector_t;
        typedef typename vector_type::scalar_type scalar_type;

        inline ray()
        {
        }
        inline ray(eClear)
        {
            clear();
        }
        inline ray(const vector_type& origin, const vector_type& direction) :
            m_origin(origin),
            m_direction(direction)
        {
        }

        inline void clear()
        {
            m_origin.clear();
            m_direction.clear();
        }

        inline const vector_type& get_origin(void) const
        {
            return m_origin;
        }
        inline void set_origin(const vector_type& origin)
        {
            m_origin = origin;
        }

        inline const vector_type& get_direction(void) const
        {
            return m_direction;
        }
        inline void set_direction(const vector_type& direction)
        {
            m_direction = direction;
        }

        inline scalar_type set_endpoints(const vector_type& start, const vector_type& end, const vector_type& def)
        {
            m_origin = start;

            m_direction = end - start;
            return m_direction.normalize(&def);
        }

        inline vector_type eval(scalar_type t) const
        {
            return m_origin + m_direction * t;
        }

    private:
        vector_type m_origin;
        vector_type m_direction;
    };

    typedef ray<vec2F> ray2F;
    typedef ray<vec3F> ray3F;

} // namespace crnlib
