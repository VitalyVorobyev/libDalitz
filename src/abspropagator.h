/** Copyright 2017 Vitaly Vorobyev
 ** @file abspropagator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_ABSPROPAGATOR_H_
#define SRC_ABSPROPAGATOR_H_

#include <complex>

class ResPropType {
 public:
    static const int NR       = 0;
    static const int RBW      = 1;
    static const int GS       = 2;
    static const int RhoOmega = 3;
    static const int Bugg     = 4;
    static const int VDst     = 5;
    static const int Flatte   = 6;
    static const int VDst2    = 7;
};

class AbsPropagator {
 public:
    AbsPropagator(const double &m, const double &p0);
    virtual ~AbsPropagator() {}
    virtual std::complex<double> operator()(const double& s,
                                            const double& p) const = 0;

    double m(void)  const {return m_m;}
    double p0(void) const {return m_p0;}

 private:
    double m_m;
    double m_p0;
};

#endif  // SRC_ABSPROPAGATOR_H_
