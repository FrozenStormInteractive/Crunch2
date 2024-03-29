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

#define CRNLIB_NO_COPY_OR_ASSIGNMENT_OP(c) \
    c(const c&); \
    c& operator=(const c&);
#define CRNLIB_NO_HEAP_ALLOC() \
private: \
    static void* operator new(size_t); \
    static void* operator new[](size_t);

namespace crnlib
{
    namespace helpers
    {
        template<typename T>
        struct rel_ops
        {
            friend bool operator!=(const T& x, const T& y)
            {
                return (!(x == y));
            }
            friend bool operator>(const T& x, const T& y)
            {
                return (y < x);
            }
            friend bool operator<=(const T& x, const T& y)
            {
                return (!(y < x));
            }
            friend bool operator>=(const T& x, const T& y)
            {
                return (!(x < y));
            }
        };

        template<typename T>
        inline T* construct(T* p)
        {
            return new (static_cast<void*>(p)) T;
        }

        template<typename T, typename U>
        inline T* construct(T* p, const U& init)
        {
            return new (static_cast<void*>(p)) T(init);
        }

        template<typename T>
        inline void construct_array(T* p, uint n)
        {
            T* q = p + n;
            for (; p != q; ++p)
            {
                new (static_cast<void*>(p)) T;
            }
        }

        template<typename T, typename U>
        inline void construct_array(T* p, uint n, const U& init)
        {
            T* q = p + n;
            for (; p != q; ++p)
            {
                new (static_cast<void*>(p)) T(init);
            }
        }

        template<typename T>
        inline void destruct(T* p)
        {
            (void)p;
            p->~T();
        }

        template<typename T>
        inline void destruct_array(T* p, uint n)
        {
            T* q = p + n;
            for (; p != q; ++p)
            {
                p->~T();
            }
        }

    } // namespace helpers
} // namespace crnlib
