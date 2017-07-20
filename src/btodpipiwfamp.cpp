#include "btodpipiwfamp.h"

#include <utility>
#include <cmath>

using compld = std::complex<double>;
using vectcd = std::vector<compld>;

using std::arg;  // complex phase
using std::move;

using std::exp;

constexpr double mB0 = 5.279;    // B0 mass
constexpr double mpi = 0.139568; // pi+ mass
constexpr double mD0 = 1.865;  // D0 mass

constexpr double fpi = 133.;  // pion decay constant
constexpr double fD = 207.;  // D0 decay constant
constexpr double FBD = 0.33;  // B -> D transition formfactor
constexpr double FBpi = 0.33;  // B -> pi transition formfactor

constexpr compld j2(0, 1);

BtoDpipiWFAmp::BtoDpipiWFAmp(double phiD0, double phiDst, double phiDstst) :
    AbsDalitzModel(mB0, mD0, mpi, mpi),
    KuzBtoDpipiAmp() {
    SetModelTitle("B0 -> D0 pi+ pi- (suppressed) decay amplitude");
    init(move(phiD0), move(phiDst), move(phiDstst));
}

// Vector of Res_amp (rho, omega and rho' combined)
void BtoDpipiWFAmp::GetResVals(vectcd* resv,
                double mABsq, double mACsq) const {
    return KuzBtoDpipiAmp::GetResVals(resv, mACsq, mABsq);
}

// Complex amplitude of a specific resonance
compld BtoDpipiWFAmp::GetResVal(double mABsq, double mACsq, int resnum) const {
    return KuzBtoDpipiAmp::GetResVal(mACsq, mABsq, resnum);
}

void BtoDpipiWFAmp::init(double phiD0, double phiDst, double phiDstst) {
    auto coefs = GetCoefficients();
    auto rhoArg = arg(coefs[3]);
    const auto coef = (fD * FBpi) / (fpi * FBD);
    coefs[0] *= coef * exp(j2 * (phiDst - rhoArg));  // D*v
    coefs[1] *= 0.1 * coef * exp(j2 * (phiDstst - rhoArg));  // D2*
    coefs[2] *= coef * exp(j2 * (phiD0 - rhoArg));  // D0*
    for (auto idx = 3; idx < 8; idx++)  // all R -> pi+ pi-
        coefs[idx] *= exp(-j2 * rhoArg);
    SetCoefficients(coefs);
}
