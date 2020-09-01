// File: ryg_dxt.hpp
#pragma once

#include "crn_ryg_types.hpp"
#include "crn_export.h"

namespace ryg_dxt {
	CRN_EXPORT extern sU8 Expand5[32];
	CRN_EXPORT extern sU8 Expand6[64];
	CRN_EXPORT extern sU8 OMatch5[256][2];
	CRN_EXPORT extern sU8 OMatch6[256][2];
	CRN_EXPORT extern sU8 OMatch5_3[256][2];
	CRN_EXPORT extern sU8 OMatch6_3[256][2];
	CRN_EXPORT extern sU8 QuantRBTab[256 + 16];
	CRN_EXPORT extern sU8 QuantGTab[256 + 16];

// initialize DXT codec. only needs to be called once.
	CRN_EXPORT void sInitDXT();

// input: a 4x4 pixel block, A8R8G8B8. you need to handle boundary cases
// yourself.
// alpha=sTRUE => use DXT5 (else use DXT1)
// quality: 0=fastest (no dither), 1=medium (dither)
	CRN_EXPORT void sCompressDXTBlock(sU8* dest, const sU32* src, sBool alpha, sInt quality);

	CRN_EXPORT void sCompressDXT5ABlock(sU8* dest, const sU32* src);

}  // namespace ryg_dxt
