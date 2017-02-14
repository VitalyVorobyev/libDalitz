/** Copyright 2017 Vitaly Vorobyev
 ** @file flattewidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_FLATTEWIDTH_H_
#define SRC_FLATTEWIDTH_H_

#include "./absvarwidth.h"

class FlatteWidth : public AbsVarWidth {
 public:
    explicit FlatteWidth(const double& m);

    double operator() (const double& s, const double& p) const;

 private:
    static const double m_g1;
    static const double m_g2;
    static const double m_pi;
    static const double m_pi_sq;
    static const double m_pi0;
    static const double m_pi0_sq;
    static const double m_K;
    static const double m_K_sq;
    static const double m_K0;
    static const double m_K0_sq;

    double rho_pipi(const double& s) const;
    double rho_KK(const double& s) const;
    double phsp_factor(const double& msq, const double&s) const;
};

#endif  // SRC_FLATTEWIDTH_H_
