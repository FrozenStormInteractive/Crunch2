#include <iostream>

#include "crnlib.h"

int main()
{
    crn_get_file_type_ext(crn_file_type::cCRNFileTypeCRN);
    std::cout << "crnlib " << crn_get_version() << std::endl;
    return 0;
}
