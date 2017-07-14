#include "btodpipiwfamp.h"

using compld = std::complex<double>;
using vectcd = std::vector<compld>;

constexpr double m_fm = 5.279;       // B0 mass
constexpr double m_pm = 0.139568;    // pi+ mass
constexpr double m_dstm  = 1.865;    // D0 mass

BtoDpipiWFAmp::BtoDpipiWFAmp(double phiD0, double phiDst) :
    AbsDalitzModel(m_fm, m_dstm, m_pm, m_pm),
    KuzBtoDpipiAmp() {
    SetModelTitle("B0 -> D0 pi+ pi- (suppressed) decay amplitude");
}
/** Vector of Res_amp (rho, omega and rho' combined) */
void BtoDpipiWFAmp::GetResVals(vectcd* resv,
                double mABsq, double mACsq) const {
    return KuzBtoDpipiAmp::GetResVals(resv, mACsq, mABsq);
}
/** Complex amplitude of a specific resonance */
compld BtoDpipiWFAmp::GetResVal(double mABsq, double mACsq,
                               const int resnum) const {
    return KuzBtoDpipiAmp::GetResVal(mACsq, mABsq, resnum);
}
