#include "crn_core.h"
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
