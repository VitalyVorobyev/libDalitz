/** Copyright 2017 Vitaly Vorobyev
 ** @file absvarwidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_ABSVARWIDTH_H_
#define SRC_ABSVARWIDTH_H_

class VarWType {
 public:
    static const int Const  = 0;
    static const int BW     = 1;
    static const int GS     = 2;
    static const int Flatte = 3;
    static const int Bugg   = 4;
};

class AbsVarWidth {
 public:
    AbsVarWidth(const double& G0, const double& m, const double& p0);
    virtual ~AbsVarWidth() {}
    virtual double operator()(const double& s, const double& p) const = 0;

    double G0(void) const {return m_G0;}
    double m(void)  const {return m_m;}
    double p0(void) const {return m_p0;}

 private:
    double m_G0;
    double m_m;
    double m_p0;
};

#endif  // SRC_ABSVARWIDTH_H_
