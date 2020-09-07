// File: crn_strutils.h
// See Copyright Notice and license at the end of inc/crnlib.h

#pragma once

#include "crn_export.h"

#ifdef WIN32
#define CRNLIB_PATH_SEPERATOR_CHAR '\\'
#else
#define CRNLIB_PATH_SEPERATOR_CHAR '/'
#endif

namespace crnlib 
{
	CRN_EXPORT char* crn_strdup(const char* pStr);
	CRN_EXPORT int crn_stricmp(const char* p, const char* q);

	CRN_EXPORT char* strcpy_safe(char* pDst, uint dst_len, const char* pSrc);

	CRN_EXPORT bool int_to_string(int value, char* pDst, uint len);
	CRN_EXPORT bool uint_to_string(uint value, char* pDst, uint len);

	CRN_EXPORT bool string_to_int(const char*& pBuf, int& value);

	CRN_EXPORT bool string_to_uint(const char*& pBuf, uint& value);

	CRN_EXPORT bool string_to_int64(const char*& pBuf, int64& value);
	CRN_EXPORT bool string_to_uint64(const char*& pBuf, uint64& value);

	CRN_EXPORT bool string_to_bool(const char* p, bool& value);

	CRN_EXPORT bool string_to_float(const char*& p, float& value, uint round_digit = 512U);

	CRN_EXPORT bool string_to_double(const char*& p, double& value, uint round_digit = 512U);
	CRN_EXPORT bool string_to_double(const char*& p, const char* pEnd, double& value, uint round_digit = 512U);
} // namespace crnlib
