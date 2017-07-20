#ifndef BTODPIPIWFAMP_H
#define BTODPIPIWFAMP_H

#include "kuzbtodpipiamp.h"

/**
 * @brief Amplitude for the suppressed B0 -> D0 pi+ pi- decay with
 * "wrong" flavor of D0 meson
 */
class BtoDpipiWFAmp : public KuzBtoDpipiAmp {
    /** Set up amplitudes */
    void init(double phiD0, double phiDst, double phiDstst);

 public:
    /**
     * @brief BtoDpipiWFAmp constructor
     * @param phiD0. Phase of the scalar D0 amplitude
     * @param phiDst. Phase of the virtual D* amplitude
     * @param phiDstst. Phase of the D*2 amplitude
     * All phases are set relative to the rho(770) resonance
     */
    BtoDpipiWFAmp(double phiD0, double phiDst, double phiDstst);

    /** Vector of Res_amp (rho, omega and rho' combined) */
    void GetResVals(std::vector<std::complex<double>>* resv,
                    double mABsq, double mACsq) const override;
    /** Complex amplitude of a specific resonance */
    std::complex<double> GetResVal(double mABsq, double mACsq,
                                   int resnum) const override;
};

#endif  // BTODPIPIWFAMP_H
