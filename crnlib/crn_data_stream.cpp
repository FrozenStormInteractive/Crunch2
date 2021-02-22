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
#include "crn_data_stream.h"

namespace crnlib
{
    data_stream::data_stream() :
        m_attribs(0),
        m_opened(false),
        m_error(false),
        m_got_cr(false)
    {
    }

    data_stream::data_stream(const char* pName, uint attribs) :
        m_name(pName),
        m_attribs(static_cast<uint16>(attribs)),
        m_opened(false),
        m_error(false),
        m_got_cr(false)
    {
    }

    uint64 data_stream::skip(uint64 len)
    {
        uint64 total_bytes_read = 0;

        const uint cBufSize = 1024;
        uint8 buf[cBufSize];

        while (len)
        {
            const uint64 bytes_to_read = math::minimum<uint64>(sizeof(buf), len);
            const uint64 bytes_read = read(buf, static_cast<uint>(bytes_to_read));
            total_bytes_read += bytes_read;

            if (bytes_read != bytes_to_read)
            {
                break;
            }

            len -= bytes_read;
        }

        return total_bytes_read;
    }

    bool data_stream::read_line(dynamic_string& str)
    {
        str.empty();

        for (;;)
        {
            const int c = read_byte();

            const bool prev_got_cr = m_got_cr;
            m_got_cr = false;

            if (c < 0)
            {
                if (!str.is_empty())
                {
                    break;
                }
                return false;
            }
            else if ((26 == c) || (!c))
            {
                continue;
            }
            else if (13 == c)
            {
                m_got_cr = true;
                break;
            }
            else if (10 == c)
            {
                if (prev_got_cr)
                {
                    continue;
                }
                break;
            }

            str.append_char(static_cast<char>(c));
        }

        return true;
    }

    bool data_stream::printf(const char* p, ...)
    {
        va_list args;

        va_start(args, p);
        dynamic_string buf;
        buf.format_args(p, args);
        va_end(args);

        return write(buf.get_ptr(), buf.get_len() * sizeof(char)) == buf.get_len() * sizeof(char);
    }

    bool data_stream::write_line(const dynamic_string& str)
    {
        if (!str.is_empty())
        {
            return write(str.get_ptr(), str.get_len()) == str.get_len();
        }

        return true;
    }

    bool data_stream::read_array(vector<uint8>& buf)
    {
        if (buf.size() < get_remaining())
        {
            if (get_remaining() > 1024U * 1024U * 1024U)
            {
                return false;
            }

            buf.resize((uint)get_remaining());
        }

        if (!get_remaining())
        {
            buf.resize(0);
            return true;
        }

        return read(&buf[0], buf.size()) == buf.size();
    }

    bool data_stream::write_array(const vector<uint8>& buf)
    {
        if (!buf.empty())
        {
            return write(&buf[0], buf.size()) == buf.size();
        }
        return true;
    }

} // namespace crnlib
