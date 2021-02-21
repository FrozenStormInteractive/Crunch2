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
    class CRN_EXPORT find_files
    {
    public:
        struct file_desc
        {
            inline file_desc() :
                m_is_dir(false)
            {
            }

            dynamic_string m_fullname;
            dynamic_string m_base;
            dynamic_string m_rel;
            dynamic_string m_name;
            bool m_is_dir;

            inline bool operator==(const file_desc& other) const
            {
                return m_fullname == other.m_fullname;
            }
            inline bool operator<(const file_desc& other) const
            {
                return m_fullname < other.m_fullname;
            }

            inline operator size_t() const
            {
                return static_cast<size_t>(m_fullname);
            }
        };

        typedef crnlib::vector<file_desc> file_desc_vec;

        inline find_files()
        {
            m_last_error = 0; // S_OK;
        }

        enum flags
        {
            cFlagRecursive = 1,
            cFlagAllowDirs = 2,
            cFlagAllowFiles = 4,
            cFlagAllowHidden = 8
        };

        bool find(const char* pBasepath, const char* pFilespec, uint flags = cFlagAllowFiles);

        bool find(const char* pSpec, uint flags = cFlagAllowFiles);

        // An HRESULT under Win32. FIXME: Abstract this better?
        inline int64 get_last_error() const
        {
            return m_last_error;
        }

        const file_desc_vec& get_files() const
        {
            return m_files;
        }

    private:
        file_desc_vec m_files;

        // A HRESULT under Win32
        int64 m_last_error;

        bool find_internal(const char* pBasepath, const char* pRelpath, const char* pFilespec, uint flags, int level);

    }; // class find_files
} // namespace crnlib
