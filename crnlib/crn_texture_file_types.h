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

#include "crnlib.h"
#include "crn_vec.h"
#include "crn_pixel_format.h"
#include "crn_export.h"

namespace crnlib
{
    struct texture_file_types
    {
        enum format
        {
            cFormatInvalid = -1,

            cFormatDDS,
            cFormatCRN,
            cFormatKTX,

            cNumMipmappedFileFormats,

            cFormatTGA = cNumMipmappedFileFormats,
            cFormatPNG,
            cFormatJPG,
            cFormatJPEG,
            cFormatBMP,
            cFormatGIF,
            cFormatTIF,
            cFormatTIFF,
            cFormatPPM,
            cFormatPGM,
            cFormatPSD,
            cFormatJP2,

            cNumRegularFileFormats,

            cNumImageFileFormats = cNumRegularFileFormats - cNumMipmappedFileFormats,

            // Not really a file format
            cFormatClipboard = cNumRegularFileFormats,
            cFormatDragDrop,

            cNumFileFormats,
        };

        CRN_EXPORT static const char* get_extension(format fmt);

        CRN_EXPORT static format determine_file_format(const char* pFilename);

        CRN_EXPORT static bool supports_mipmaps(format fmt);
        CRN_EXPORT static bool supports_alpha(format fmt);
    };

    enum texture_type
    {
        cTextureTypeUnknown = 0,
        cTextureTypeRegularMap,
        cTextureTypeNormalMap,
        cTextureTypeVerticalCrossCubemap,
        cTextureTypeCubemap,

        cNumTextureTypes
    };

    CRN_EXPORT const char* get_texture_type_desc(texture_type t);
} // namespace crnlib
