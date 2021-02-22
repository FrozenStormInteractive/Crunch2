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
#include "crn_texture_file_types.h"
#include "crn_file_utils.h"

namespace crnlib
{
    const char* texture_file_types::get_extension(format fmt)
    {
        CRNLIB_ASSERT(fmt < cNumFileFormats);
        if (fmt >= cNumFileFormats)
        {
            return nullptr;
        }

        static const char* extensions[cNumFileFormats] = {
            "dds",
            "crn",
            "ktx",

            "tga",
            "png",
            "jpg",
            "jpeg",
            "bmp",
            "gif",
            "tif",
            "tiff",
            "ppm",
            "pgm",
            "psd",
            "jp2",

            "<clipboard>",
            "<dragdrop>"
        };
        return extensions[fmt];
    }

    texture_file_types::format texture_file_types::determine_file_format(const char* pFilename)
    {
        dynamic_string ext;
        if (!file_utils::split_path(pFilename, nullptr, nullptr, nullptr, &ext))
        {
            return cFormatInvalid;
        }

        if (ext.is_empty())
        {
            return cFormatInvalid;
        }

        if (ext[0] == '.')
        {
            ext.right(1);
        }

        for (uint i = 0; i < cNumFileFormats; i++)
        {
            if (ext == get_extension(static_cast<format>(i)))
            {
                return static_cast<format>(i);
            }
        }
        return cFormatInvalid;
    }

    bool texture_file_types::supports_mipmaps(format fmt)
    {
        switch (fmt)
        {
        case cFormatCRN:
        case cFormatDDS:
        case cFormatKTX:
            return true;
        default:
            break;
        }

        return false;
    }

    bool texture_file_types::supports_alpha(format fmt)
    {
        switch (fmt)
        {
        case cFormatJPG:
        case cFormatJPEG:
        case cFormatGIF:
        case cFormatJP2:
            return false;
        default:
            break;
        }

        return true;
    }

    const char* get_texture_type_desc(texture_type t)
    {
        switch (t)
        {
        case cTextureTypeUnknown:
            return "Unknown";
        case cTextureTypeRegularMap:
            return "2D map";
        case cTextureTypeNormalMap:
            return "Normal map";
        case cTextureTypeVerticalCrossCubemap:
            return "Vertical Cross Cubemap";
        case cTextureTypeCubemap:
            return "Cubemap";
        default:
            break;
        }

        CRNLIB_ASSERT(false);

        return "?";
    }
} // namespace crnlib
