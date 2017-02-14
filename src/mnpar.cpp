/** Copyright 2017 Vitaly Vorobyev
 ** @file mnpar.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/mnpar.h"

using std::to_string;

typedef std::string str;

MnPar::MnPar(const str& n, const double& iv, const double& ie,
             const double& ll, const double& rl, bool fix) :
    name(n), ival(iv), ierr(ie), rlim(rl), llim(ll), fixed(fix) {}

const char* MnPar::CName(void) const {
    return name.c_str();
}

str MnPar::Parstr(void) const {
    str text = name + ": " + to_string(ival) + " +- " + to_string(ival) +
            " [" + to_string(llim) + "," + to_string(rlim) + "]";
    if (fixed) text += " (fixed)";
    return text;
}
