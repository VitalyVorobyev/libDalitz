/** Copyright 2017 Vitaly Vorobyev
 ** @file kuzbtodpipiamp.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_KUZBTODPIPIAMP_H_
#define SRC_KUZBTODPIPIAMP_H_

#include <vector>
#include <complex>
#include <string>

#include "./abssymdalitzmodel.h"

/** \brief The KuzBtoDpipiAmp class reproduceses Alex Kuzmin's
 *  Fortran code for the B0 -> D0 pi+ pi- decay amplitude model
 *  describing the Belle experiment data. Ref!!! **/
class KuzBtoDpipiAmp : public AbsSymDalitzModel {
    void init(void);

 public:
    KuzBtoDpipiAmp();
    /** Vector of Res_amp (rho, omega and rho' combined) */
    void GetResVals(std::vector<std::complex<double>>* resv,
                    double mABsq, double mACsq) const override;
    /** Complex amplitude of a specific resonance */
    std::complex<double> GetResVal(double mABsq, double mACsq,
                                   const int resnum) const override;
    /** Set parameters */
    void SetParams(const std::vector<double>& pars) override;
    /** Get normalization */
    double Norm(void) const override;
    /** Print info in standard output */
    void PrintSummary(void);
};

#endif  // SRC_KUZBTODPIPIAMP_H_
