#ifndef ABSVARWIDTH_H
#define ABSVARWIDTH_H

class VarWType{
public:
  static const int Const = 0;
  static const int BW = 1;
  static const int GS = 2;
  static const int Flatte = 3;
  static const int Bugg = 4;
};

class AbsVarWidth{
public:
  AbsVarWidth(const double& G0, const double& m, const double& p0);
  virtual double operator()(const double& s, const double& p) const = 0;

  double G0(void) const {return m_G0;}
  double m(void)  const {return m_m;}
  double p0(void) const {return m_p0;}

private:
  double m_G0;
  double m_m;
  double m_p0;
};

#endif // ABSVARWIDTH_H
