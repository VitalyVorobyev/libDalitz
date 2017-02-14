/** Copyright 2017 Vitaly Vorobyev
 ** @file mnpar.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_MNPAR_H_
#define SRC_MNPAR_H_

#include <string>

///
/// \brief The MnPar class
///
class MnPar {
 public:
    MnPar(const std::string& n, const double& iv, const double& ie,
          const double& ll, const double& rl, bool fix = true);

    void Release(void) {fixed = false;}
    void Fix(void) {fixed = true;}

    const char* CName(void) const;
    std::string Parstr(void) const;

    std::string name;
    double ival;
    double ierr;
    double rlim;
    double llim;
    bool fixed;
};

#endif  // SRC_MNPAR_H_
