/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzmcintegral.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_DALITZMCINTEGRAL_H_
#define SRC_DALITZMCINTEGRAL_H_

#include <string>
#include <vector>
#include <complex>
#include <cstdint>

#include "./randomdalitzpoint.h"
#include "./absdalitzmodel.h"

///
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
///
class DalitzMCIntegral : public RandomDalitzPoint {
 public:
    explicit DalitzMCIntegral(const AbsDalitzModel* model);

    /// MC calculation of the Dalirz plot area
    double GetDalitzPlotArea(const uint64_t& nc = 0) const;
    /// Get normalization integral for DalitzModel
    int CalcIntegral(double* val, double* err, const uint64_t &nc = 0,
                     const bool morm_flag = false) const;
    /// \brief CalcBranchings calculates fit fractions of each
    /// resonance in Dalitz model
    double CalcBranchings(std::vector<double>* brvec,
                          std::vector<double>* brerr,
                          const uint64_t& nc = 0) const;
    /// Calculates all integrals needed for fast calculation of normalization
    int CalcNormMap(std::string* ofile, const uint64_t& nc = 0) const;
    /// \brief DalitzMCIntegral::SmartNromMap
    std::vector<std::vector<std::complex<double>>>
    SmartNormMap(const uint64_t& nc = 0) const;
    ///
    void CheckNormalization(const uint64_t &nc = 0) const;
    /// \brief NCounts
    uint64_t NCounts(void) const {return m_ncounts;}
    /// \brief SetNCounts
    void SetNCounts(const uint64_t& p) {m_ncounts = p; return;}
    ///
    void SetPrefix(const std::string& prefix) {m_prefix = prefix;}

 private:
    double GetGaussIntegral(const unsigned resnum, const uint64_t &nc) const;
    std::complex<double> GetGaussIntegral(const unsigned resnum1,
                                          const unsigned resnum2,
                                          const uint64_t& nc) const;
    void ShowIntMatrix(
            const std::vector<std::vector<std::complex<double>>>& vals) const;
    void ShowIntMatrix(
            const std::vector<std::vector<std::complex<double>>>& vals,
            const std::vector<std::vector<std::complex<double>>>& errs) const;
    int WriteIntMatrix(
            const std::string& fname,
            const std::vector<std::vector<std::complex<double>>>& vals) const;
    int WriteIntMatrix(
            const std::string& fname,
            const std::vector<std::vector<std::complex<double>>>& vals,
            const std::vector<std::vector<std::complex<double>>>& errs) const;

    const AbsDalitzModel* m_model;
    uint64_t m_ncounts;
    std::string m_prefix;

    static double normal_pdf(const double x, const double m, const double s);
    static const double inv_sqrt_2pi;
};

#endif  // SRC_DALITZMCINTEGRAL_H_

