// File: crn_checksum.h
#pragma once

#include "crn_export.h"

namespace crnlib {
const uint cInitAdler32 = 1U;
CRN_EXPORT uint adler32(const void* pBuf, size_t buflen, uint adler32 = cInitAdler32);

// crc16() intended for small buffers - doesn't use an acceleration table.
const uint cInitCRC16 = 0;
CRN_EXPORT uint16 crc16(const void* pBuf, size_t len, uint16 crc = cInitCRC16);

}  // namespace crnlib
