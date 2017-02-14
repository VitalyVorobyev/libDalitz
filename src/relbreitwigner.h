/** Copyright 2017 Vitaly Vorobyev
 ** @file relbreitwigner.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_RELBREITWIGNER_H_
#define SRC_RELBREITWIGNER_H_

#include <complex>

#include "./abspropagator.h"
#include "./absvarwidth.h"

///
/// \brief Class implementing relativistic Breit-Wigner lineshape
/// There are three options for variable resonance width:
///  1. Constant width
///  2. Standard for RBW dependence
///  3. Flatte
class RelBreitWigner : public AbsPropagator {
 public:
    RelBreitWigner(const double& G0, const double& m, const double& p0,
                   const int mom, const int wtype = VarWType::BW);
    ~RelBreitWigner();

    std::complex<double> operator()(const double& s, const double& p) const;

 private:
    int m_wtype;
    AbsVarWidth* m_width;
};

#endif  // SRC_RELBREITWIGNER_H_
