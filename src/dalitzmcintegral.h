#ifndef DALITZMCINTEGRAL_H
#define DALITZMCINTEGRAL_H

#include "dalitzmodel.h"
#include "randomdalitzpoint.h"

/// \brief Class for MC integration of DalitzModel.
/// Method CalcIntegral() calculates normalization integral usefull in
/// maximum likelihood fit.
///
/// Method CalcBranchings() calculates branching
/// fractions of each DalitzPlotObject in DalitzModel.
///
/// Integration method is based on formulas at
/// http://mathworld.wolfram.com/MonteCarloIntegration.html
///
/// int f dV = V<f> + V*sqrt((<f^2>-<f>^2)/N), where
/// <f>   = (sum_i f(x_i))/N and
/// <f^2> = (sum_i f^2(x_i))/N and
/// V is square of Dalitz plot.
///
/// See Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling,
/// W. T. "Simple Monte Carlo Integration" and "Adaptive and Recursive Monte
/// Carlo Methods." §7.6 and 7.8 in Numerical Recipes in FORTRAN:
/// The Art of Scientific Computing, 2nd ed. Cambridge, England: Cambridge
/// University Press, pp. 295-299 and 306-319, 1992.
///
/// See also Stefan Weinzierl, "Introduction to Monte Carlo methods",
/// http://arxiv.org/abs/hep-ph/0006269

class DalitzMCIntegral : public RandomDalitzPoint{
public:
  DalitzMCIntegral(const DalitzModel& model);
  int CalcIntegral(double& val, double& err, const long& nc = 0); /// Get normalization integral for DalitzModel
  double CalcBranchings(std::vector<double>& brvec,std::vector<double>& brerr,const long& nc = 0); /// Get branching fractions of each DalitzPlotObject in DalitzModel

  void SetNCounts(const long& p){m_ncounts = p; return;}
  long GetNCounts(void) const {return m_ncounts;}
private:
  DalitzModel* m_model;
  long m_ncounts;
  double m_int;
  double m_err;
};

#endif // DALITZMCINTEGRAL_H
