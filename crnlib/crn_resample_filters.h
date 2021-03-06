// RG: This is public domain code, originally derived from Graphics Gems 3, see: http://code.google.com/p/imageresampler/

#pragma once

#include "crn_export.h"

namespace crnlib
{
    typedef float (*resample_filter_func)(float t);

    struct resample_filter
    {
        char name[32];
        resample_filter_func func;
        float support;
    };

    CRN_EXPORT extern const resample_filter g_resample_filters[];
    CRN_EXPORT extern const int g_num_resample_filters;

    CRN_EXPORT int find_resample_filter(const char* pName);
} // namespace crnlib
