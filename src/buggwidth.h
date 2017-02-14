/** Copyright 2017 Vitaly Vorobyev
 ** @file buggwidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BUGGWIDTH_H_
#define SRC_BUGGWIDTH_H_

#include "./absvarwidth.h"

/**
 * @brief The BuggWidth class. D.V. Bugg, "The mass of the σ pole",
 * Journal of Physics G: Nuclear and Particle Physics,
 * vol 34, num 1, page 151, year 2007
 * http://stacks.iop.org/0954-3899/34/i=1/a=011
 */
class BuggWidth : public AbsVarWidth {
 public:
    BuggWidth(void);
    double mrGamma1(const double& s);
    void GetWidths(const double& s, double* G1, double* GTot);
    double sA(void) const {return m_sA;}
    double mrsq(void) const {return m_mrsq;}
    double g1sq(void) const {return m_g1sq_pc;}
    double z(void) const {return m_z_pc;}

    double operator()(const double& s = 0, const double& p = 0) const;

 private:
    static const double m_mr;
    static const double m_mrsq;
    static const double m_sA;
    static const double m_b1;
    static const double m_b2;
    static const double m_A;
    static const double m_g4pi;
    static const double m_alpha;
    static const double m_PI_Mass;
    static const double m_K_Mass;
    static const double m_eta_Mass;

    double m_g1sq_pc;
    double m_rho1_ps;
    double m_z_pc;

    inline double rho_4pi(const double& s);
    inline double rho(const double& m, const double& s);

    double j1(const double& s);
    inline double z(const double& s);
    inline double g1sq(const double& s);

    double mrGamma2(const double& s);
    double mrGamma3(const double& s);
    double mrGamma4(const double& s);
};

#endif  // SRC_BUGGWIDTH_H_
