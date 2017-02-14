/** Copyright 2017 Vitaly Vorobyev
 ** @file formfactor.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_FORMFACTOR_H_
#define SRC_FORMFACTOR_H_

class FormFactor {
 public:
    FormFactor(const double &_r, const double &_p0);
    virtual double operator() (const double& p) const = 0;

    double r(void)  const {return m_r;}
    double p0(void) const {return m_p0;}

 private:
    double m_r;  // resonance radial parameter
    double m_p0;
};

#endif  // SRC_FORMFACTOR_H_
