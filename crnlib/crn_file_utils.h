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
    struct CRN_EXPORT file_utils
    {
        // Returns true if pSrcFilename is older than pDstFilename
        static bool is_read_only(const char* pFilename);
        static bool disable_read_only(const char* pFilename);
        static bool is_older_than(const char* pSrcFilename, const char* pDstFilename);
        static bool does_file_exist(const char* pFilename);
        static bool does_dir_exist(const char* pDir);
        static bool get_file_size(const char* pFilename, uint64& file_size);
        static bool get_file_size(const char* pFilename, uint32& file_size);

        static bool is_path_separator(char c);
        static bool is_path_or_drive_separator(char c);
        static bool is_drive_separator(char c);

        static bool split_path(const char* p, dynamic_string* pDrive, dynamic_string* pDir, dynamic_string* pFilename, dynamic_string* pExt);
        static bool split_path(const char* p, dynamic_string& path, dynamic_string& filename);

        static bool get_pathname(const char* p, dynamic_string& path);
        static bool get_filename(const char* p, dynamic_string& filename);

        static void combine_path(dynamic_string& dst, const char* pA, const char* pB);
        static void combine_path(dynamic_string& dst, const char* pA, const char* pB, const char* pC);

        static bool full_path(dynamic_string& path);
        static bool get_extension(dynamic_string& filename);
        static bool remove_extension(dynamic_string& filename);
        static bool create_path(const dynamic_string& path);
        static void trim_trailing_seperator(dynamic_string& path);

        static int wildcmp(const char* pWild, const char* pString);

        static bool write_buf_to_file(const char* pPath, const void* pData, size_t data_size);

    }; // struct file_utils
} // namespace crnlib
