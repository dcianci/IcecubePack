#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <getopt.h>

#include "nuSQuIDS/nuSQuIDS.h"

std::array < std::array < std::array < float,200 >,21 >, 10 > syst_resp_array_nu_high, syst_resp_array_nubar_high, syst_resp_array_nu_low, syst_resp_array_nubar_low;

std::array < std::array < std::array < float,200 >,21 >, 10 > interpSystRespArray_nu(float eff);
std::array < std::array < std::array < float,200 >,21 >, 10 > interpSystRespArray_nubar(float eff);
