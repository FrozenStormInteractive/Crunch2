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

#ifndef CRNLIB_MIN_ALLOC_ALIGNMENT
#define CRNLIB_MIN_ALLOC_ALIGNMENT sizeof(size_t) * 2
#endif

namespace crnlib
{
#if CRNLIB_64BIT_POINTERS
    const uint64 CRNLIB_MAX_POSSIBLE_BLOCK_SIZE = 0x400000000ULL;
#else
    const uint32 CRNLIB_MAX_POSSIBLE_BLOCK_SIZE = 0x7FFF0000U;
#endif

    CRN_EXPORT void* crnlib_malloc(size_t size);
    CRN_EXPORT void* crnlib_malloc(size_t size, size_t* pActual_size);
    CRN_EXPORT void* crnlib_realloc(void* p, size_t size, size_t* pActual_size = nullptr, bool movable = true);
    CRN_EXPORT void* crnlib_calloc(size_t count, size_t size, size_t* pActual_size = nullptr);
    CRN_EXPORT void crnlib_free(void* p);
    CRN_EXPORT size_t crnlib_msize(void* p);
    CRN_EXPORT void crnlib_print_mem_stats();
    CRN_EXPORT void crnlib_mem_error(const char* p_msg);

    // omfg - there must be a better way

    template<typename T>
    inline T* crnlib_new()
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        if (CRNLIB_IS_SCALAR_TYPE(T))
        {
            return p;
        }
        return helpers::construct(p);
    }

    template<typename T, typename A>
    inline T* crnlib_new(const A& init0)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0);
    }

    template<typename T, typename A>
    inline T* crnlib_new(A& init0)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0);
    }

    template<typename T, typename A, typename B>
    inline T* crnlib_new(const A& init0, const B& init1)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1);
    }

    template<typename T, typename A, typename B, typename C>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2);
    }

    template<typename T, typename A, typename B, typename C, typename D>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6, const H& init7)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6, init7);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6, const H& init7, const I& init8)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6, init7, init8);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6, const H& init7, const I& init8, const J& init9)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6, init7, init8, init9);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J, typename K>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6, const H& init7, const I& init8, const J& init9, const K& init10)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6, init7, init8, init9, init10);
    }

    template<typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J, typename K, typename L>
    inline T* crnlib_new(const A& init0, const B& init1, const C& init2, const D& init3, const E& init4, const F& init5, const G& init6, const H& init7, const I& init8, const J& init9, const K& init10, const L& init11)
    {
        T* p = static_cast<T*>(crnlib_malloc(sizeof(T)));
        return new (static_cast<void*>(p)) T(init0, init1, init2, init3, init4, init5, init6, init7, init8, init9, init10, init11);
    }

    template<typename T>
    inline T* crnlib_new_array(uint32 num)
    {
        if (!num)
        {
            num = 1;
        }

        uint64 total = CRNLIB_MIN_ALLOC_ALIGNMENT + sizeof(T) * num;
        if (total > CRNLIB_MAX_POSSIBLE_BLOCK_SIZE)
        {
            crnlib_mem_error("crnlib_new_array: Array too large!");
            return nullptr;
        }
        uint8* q = static_cast<uint8*>(crnlib_malloc(static_cast<size_t>(total)));

        T* p = reinterpret_cast<T*>(q + CRNLIB_MIN_ALLOC_ALIGNMENT);

        reinterpret_cast<uint32*>(p)[-1] = num;
        reinterpret_cast<uint32*>(p)[-2] = ~num;

        if (!CRNLIB_IS_SCALAR_TYPE(T))
        {
            helpers::construct_array(p, num);
        }
        return p;
    }

    template<typename T>
    inline void crnlib_delete(T* p)
    {
        if (p)
        {
            if (!CRNLIB_IS_SCALAR_TYPE(T))
            {
                helpers::destruct(p);
            }
            crnlib_free(p);
        }
    }

    template<typename T>
    inline void crnlib_delete_array(T* p)
    {
        if (p)
        {
            const uint32 num = reinterpret_cast<uint32*>(p)[-1];
            const uint32 num_check = reinterpret_cast<uint32*>(p)[-2];
            CRNLIB_ASSERT(num && (num == ~num_check));
            if (num == ~num_check)
            {
                if (!CRNLIB_IS_SCALAR_TYPE(T))
                {
                    helpers::destruct_array(p, num);
                }

                crnlib_free(reinterpret_cast<uint8*>(p) - CRNLIB_MIN_ALLOC_ALIGNMENT);
            }
        }
    }

} // namespace crnlib

#define CRNLIB_DEFINE_NEW_DELETE \
    void* operator new(size_t size) \
    { \
        void* p = crnlib::crnlib_malloc(size); \
        if (!p) \
        { \
            crnlib_fail("new: Out of memory!", __FILE__, __LINE__); \
        } \
        return p; \
    } \
    void* operator new[](size_t size) \
    { \
        void* p = crnlib::crnlib_malloc(size); \
        if (!p) \
        { \
            crnlib_fail("new[]: Out of memory!", __FILE__, __LINE__); \
        } \
        return p; \
    } \
    void operator delete(void* p_block) \
    { \
        crnlib::crnlib_free(p_block); \
    } \
    void operator delete[](void* p_block) \
    { \
        crnlib::crnlib_free(p_block); \
    }
