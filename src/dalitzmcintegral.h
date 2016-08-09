#ifndef DALITZMCINTEGRAL_H
#define DALITZMCINTEGRAL_H

#include "randomdalitzpoint.h"
#include "absdalitzmodel.h"

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
/// V is area of Dalitz plot.
///
/// See Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling,
/// W. T. "Simple Monte Carlo Integration" and "Adaptive and Recursive Monte
/// Carlo Methods." §7.6 and 7.8 in Numerical Recipes in FORTRAN:
/// The Art of Scientific Computing, 2nd ed. Cambridge, England: Cambridge
/// University Press, pp. 295-299 and 306-319, 1992.
///
/// See also Stefan Weinzierl, "Introduction to Monte Carlo methods",
/// http://arxiv.org/abs/hep-ph/0006269

typedef unsigned long long U64;

class DalitzMCIntegral : public RandomDalitzPoint{
public:
  DalitzMCIntegral(const AbsDalitzModel* model);

  /// MC calculation of the Dalirz plot area
  double GetDalitzPlotArea(const U64& nc = 0) const;
  /// Get normalization integral for DalitzModel
  int CalcIntegral(double& val, double& err, const U64 &nc = 0, const bool morm_flag = false) const;
  /// \brief CalcBranchings calculates fit fractions of each resonance in Dalitz model
  double CalcBranchings(vectd& brvec,vectd& brerr,const U64& nc = 0) const;
  /// Calculates all integrals needed for fast calculation of normalization
  std::vector<vectcd> CalcNormMap(str &ofile, const U64& nc = 0) const;
  /// \brief DalitzMCIntegral::SmartNromMap
  std::vector<vectcd> SmartNormMap(const U64& nc = 0) const;
  ///
  void CheckNormalization(const U64 &nc = 0) const;
  /// \brief NCounts
  long NCounts(void) const {return m_ncounts;}
  /// \brief SetNCounts
  void SetNCounts(const long& p){m_ncounts = p; return;}
  ///
  void SetPrefix(cstr& prefix) {m_prefix = prefix;}
private:
  double GetGaussIntegral(const int resnum, const U64 &nc) const;
  compld GetGaussIntegral(const int resnum1,const int resnum2, const U64& nc) const;
  void ShowIntMatrix(const std::vector<vectcd>& vals) const;
  void ShowIntMatrix(const std::vector<vectcd>& vals,const std::vector<vectcd>& errs) const;
  int WriteIntMatrix(const str& fname,const std::vector<vectcd>& vals) const;
  int WriteIntMatrix(const str& fname,const std::vector<vectcd>& vals,const std::vector<vectcd>& errs) const;

  const AbsDalitzModel* m_model;
  U64 m_ncounts;
  str m_prefix;
};

#endif // DALITZMCINTEGRAL_H

