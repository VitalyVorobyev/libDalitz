#ifndef BTODPIPIWFAMP_H
#define BTODPIPIWFAMP_H

#include "kuzbtodpipiamp.h"

/**
 * @brief Amplitude for the suppressed B0 -> D0 pi+ pi- decay with
 * "wrong" flavor of D0 meson
 */
class BtoDpipiWFAmp : public KuzBtoDpipiAmp {
public:
    BtoDpipiWFAmp(double phiD0, double phiDst);

    /** Vector of Res_amp (rho, omega and rho' combined) */
    void GetResVals(std::vector<std::complex<double>>* resv,
                    double mABsq, double mACsq) const override;
    /** Complex amplitude of a specific resonance */
    std::complex<double> GetResVal(double mABsq, double mACsq,
                                   const int resnum) const override;
};

#endif  // BTODPIPIWFAMP_H
