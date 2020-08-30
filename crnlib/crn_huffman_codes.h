// File: crn_huffman_codes.h
// See Copyright Notice and license at the end of inc/crnlib.h
#pragma once

#include "crn_export.h"

namespace crnlib {
const uint cHuffmanMaxSupportedSyms = 8192;

CRN_EXPORT void* create_generate_huffman_codes_tables();
CRN_EXPORT void free_generate_huffman_codes_tables(void* p);

CRN_EXPORT bool generate_huffman_codes(void* pContext, uint num_syms, const uint16* pFreq, uint8* pCodesizes, uint& max_code_size, uint& total_freq_ret);

}  // namespace crnlib
