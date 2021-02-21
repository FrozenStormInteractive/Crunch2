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
    struct Area
    {
        struct Area *Pprev, *Pnext;

        int x1, y1, x2, y2;

        uint get_width() const
        {
            return x2 - x1 + 1;
        }
        uint get_height() const
        {
            return y2 - y1 + 1;
        }
        uint get_area() const
        {
            return get_width() * get_height();
        }
    };

    typedef Area* Area_Ptr;

    struct Area_List
    {
        int total_areas;
        int next_free;

        Area *Phead, *Ptail, *Pfree;
    };

    typedef Area_List* Area_List_Ptr;

    CRN_EXPORT Area_List* Area_List_init(int max_areas);
    CRN_EXPORT void Area_List_deinit(Area_List* Pobj_base);

    CRN_EXPORT void Area_List_print(Area_List* Plist);

    CRN_EXPORT Area_List* Area_List_dup_new(Area_List* Plist, int x_ofs, int y_ofs);

    CRN_EXPORT uint Area_List_get_num(Area_List* Plist);

    // src and dst area lists must have the same number of total areas.
    CRN_EXPORT void Area_List_dup(Area_List* Psrc_list, Area_List* Pdst_list, int x_ofs, int y_ofs);

    CRN_EXPORT void Area_List_copy(Area_List* Psrc_list, Area_List* Pdst_list, int x_ofs, int y_ofs);

    CRN_EXPORT void Area_List_clear(Area_List* Plist);

    CRN_EXPORT void Area_List_set(Area_List* Plist, int x1, int y1, int x2, int y2);

    // logical: x and (not y)
    CRN_EXPORT void Area_List_remove(Area_List* Plist, int x1, int y1, int x2, int y2);

    // logical: x or y
    CRN_EXPORT void Area_List_insert(Area_List* Plist, int x1, int y1, int x2, int y2, bool combine);

    // logical: x and y
    CRN_EXPORT void Area_List_intersect_area(Area_List* Plist, int x1, int y1, int x2, int y2);

    // logical: x and y
    CRN_EXPORT void Area_List_intersect_Area_List(Area_List* Pouter_list, Area_List* Pinner_list, Area_List* Pdst_list);

    CRN_EXPORT Area_List_Ptr Area_List_create_optimal(Area_List_Ptr Plist);

} // namespace crnlib
