#ifndef FORMFACTOR_H
#define FORMFACTOR_H

class FormFactor{
public:
  FormFactor(const double &_r, const double &_p0sq);
  virtual double operator()(const double& p) const = 0;

  double r(void)  const {return m_r;}
  double p0sq(void) const {return m_p0sq;}

private:
  double m_r; // resonance radial parameter
  double m_p0sq;
};

#endif // FORMFACTOR_H
