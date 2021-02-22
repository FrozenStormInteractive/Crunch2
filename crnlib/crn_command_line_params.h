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

#include "crn_value.h"
#include <map>
#include "crn_export.h"

namespace crnlib
{
    // Returns the command line passed to the app as a string.
    // On systems where this isn't trivial, this function combines together the separate arguments, quoting and adding spaces as needed.
    CRN_EXPORT void get_command_line_as_single_string(dynamic_string& cmd_line, int argc, char* argv[]);

    class CRN_EXPORT command_line_params
    {
    public:
        struct param_value
        {
            inline param_value() :
                m_index(0),
                m_modifier(0)
            {
            }

            dynamic_string_array m_values;
            uint m_index;
            int8 m_modifier;
        };

        typedef std::multimap<dynamic_string, param_value> param_map;
        typedef param_map::const_iterator param_map_const_iterator;
        typedef param_map::iterator param_map_iterator;

        command_line_params();

        void clear();

        static bool split_params(const char* p, dynamic_string_array& params);

        struct param_desc
        {
            const char* m_pName;
            uint m_num_values;
            bool m_support_listing_file;
        };

        bool parse(const dynamic_string_array& params, uint n, const param_desc* pParam_desc);
        bool parse(const char* pCmd_line, uint n, const param_desc* pParam_desc, bool skip_first_param = true);

        const dynamic_string_array& get_array() const
        {
            return m_params;
        }

        bool is_param(uint index) const;

        const param_map& get_map() const
        {
            return m_param_map;
        }

        uint get_num_params() const
        {
            return static_cast<uint>(m_param_map.size());
        }

        param_map_const_iterator begin() const
        {
            return m_param_map.begin();
        }
        param_map_const_iterator end() const
        {
            return m_param_map.end();
        }

        uint find(uint num_keys, const char** ppKeys, crnlib::vector<param_map_const_iterator>* pIterators, crnlib::vector<uint>* pUnmatched_indices) const;

        void find(const char* pKey, param_map_const_iterator& begin, param_map_const_iterator& end) const;

        uint get_count(const char* pKey) const;

        // Returns end() if param cannot be found, or index is out of range.
        param_map_const_iterator get_param(const char* pKey, uint index) const;

        bool has_key(const char* pKey) const
        {
            return get_param(pKey, 0) != end();
        }

        bool has_value(const char* pKey, uint index) const;
        uint get_num_values(const char* pKey, uint index) const;

        bool get_value_as_bool(const char* pKey, uint index = 0, bool def = false) const;

        int get_value_as_int(const char* pKey, uint index, int def, int l = INT_MIN, int h = INT_MAX, uint value_index = 0) const;
        float get_value_as_float(const char* pKey, uint index, float def = 0.0f, float l = -math::cNearlyInfinite, float h = math::cNearlyInfinite, uint value_index = 0) const;

        bool get_value_as_string(const char* pKey, uint index, dynamic_string& value, uint value_index = 0) const;
        const dynamic_string& get_value_as_string_or_empty(const char* pKey, uint index = 0, uint value_index = 0) const;

    private:
        dynamic_string_array m_params;

        param_map m_param_map;

        static bool load_string_file(const char* pFilename, dynamic_string_array& strings);
    };

} // namespace crnlib
