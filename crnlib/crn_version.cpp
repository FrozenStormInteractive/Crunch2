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

#include "crnlib.h"
#include "crn_version.h"

const char* crn_get_version()
{
    return CRN_VERSION_STR;
}

int crn_get_version_number()
{
    return CRN_VERSION_NUMBER;
}

int crn_get_version_major()
{
    return CRN_VERSION_MAJOR;
}

int crn_get_version_minor()
{
    return CRN_VERSION_MINOR;
}

int crn_get_version_patch()
{
    return CRN_VERSION_PATCH;
}
