/** Copyright 2017 Vitaly Vorobyev
 ** @file kuzbtodpipiamp.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/kuzbtodpipiamp.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>

using KAmp = KuzBtoDpipiAmp;
using compld = std::complex<double>;
using vectd = std::vector<double>;
using vecti = std::vector<int>;
using vectcd = std::vector<compld>;

constexpr compld ione(0, 1);

using std::exp;
using std::pow;
using std::sqrt;

using std::cout;
using std::cerr;
using std::endl;

using std::make_pair;
using std::to_string;

/** Energy of particle A in the M -> AB decay in the M rest frame */
inline auto enrg(double M, double A2, double B2) {
    return 0.5 * (M * M + A2 - B2) / M;
}
/** A particle momentum from its energy and mass */
inline auto momt(double esq, double msq) {
    return (esq - msq) > 0 ? sqrt(esq - msq) : 0.;
}
/** Momentum of particle A in the M -> AB decay in the M rest frame */
auto momt(double M, double A2, double B2) {
    auto e = enrg(M, A2, B2);
    return momt(e * e, A2);
}

// Raw amplitudes and phases
constexpr double m_a0 = 0.35;        // par  1: D0 amplitude
constexpr double m_f0 = -2.88;       // par  3: D0 phase
constexpr double m_a1 = 0.065;       // par  2: D* amplitude
constexpr double m_f1 = -2.53;       // par  4: D* phase
constexpr double m_ab = 0;           // par 10: f(980) amplitude
constexpr double m_fb = -3.07;       // par 11: f(980) phase
constexpr double m_a3 = 1.961500;    // par 24: f(1370) amplitude
constexpr double m_f3 = -2.43;       // par 25: f(1370) phase
constexpr double m_ar = 0.5;         // par 14: rho(770) amplitude
constexpr double m_fr = 1.81;        // par 15: rho(770) phase
constexpr double m_af0 = 0.1;        // par 18: f0(550) amplitude
constexpr double m_ff0 = -0.40;      // par 19: f0(550) phase relative
                                     //         to the rho(770) phase
constexpr double m_af2 = 0.12;       // par 22: f2 amplitude
constexpr double m_ff2 = 2.77;       // par 23: f2 phase relative to
                                     //         the rho(770) phase
constexpr double m_amom = 0.028;     // par 30: omega amplitude

// Masses and widths
constexpr double m_fm = 5.279;       // B0 mass
constexpr double m_pm = 0.139568;    // pi+ mass
constexpr double m_dstm  = 1.865;    // D0 mass
constexpr double m_amf0 = 0.51300;   // par 16: f0 mass
constexpr double m_gmf0 = 0.33500;   // par 17: f0 width
constexpr double m_amf2 = 1.2750;    // par 20: f2 mass
constexpr double m_gmf2 = 0.18500;   // par 21: f2 width
constexpr double m_amf098 = 0.97;    // f(980) mass
constexpr double m_gmf098 = 0.0440;  // f(980) width
constexpr double m_am3 = 1.434;      // f(1370) mass
constexpr double m_gm3 = 0.173;      // f(1370) width
constexpr double m_ams0 = 2.3080;    // par 7: D*0 mass
constexpr double m_gms0 = 0.27611;   // par 8: D*0 width
constexpr double m_ams2 = 2.4677;    // par 5: D2 mass
constexpr double m_gms2 = 0.056002;  // par 6: D2 width
constexpr double m_amrh = 0.77560;   // par 12: rho(770) mass
constexpr double m_gmrh = 0.15000;   // par 13: rho(770) width
constexpr double m_amrh1 = 1.465;    // Mass rho'
constexpr double m_gmrh1 = 0.31;     // Width rho'

constexpr double m_rhop_amp = -0.248; // rho' amplitude relative to rho(770)

// Auxiliary parameters
constexpr double m_fm2 = m_fm * m_fm;  // B0 mass sq
constexpr double m_pm2   = m_pm * m_pm;  // pi+ mass sq
constexpr double m_dstm2 = m_dstm * m_dstm;  // D0 mass sq
constexpr double m_ddd = 1. - m_a0 + m_a1 + m_ab + m_a3 + m_ar + m_af0 + m_af2;

inline double dot(const vectd& u, const vectd& v) {
    return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
}

// Amplitude wrappers
auto anf0(double am, double ag) {
    static vectd pars = {-69.398, 608.83, -1152.1, 940.13,
                         -198.62, -148.03, 66.298};
    vectd vars = {1.};
    for (auto idx = 1u; idx < pars.size(); idx++)
        vars.emplace_back(vars[idx-1] * am);
    return dot(pars, vars) / (pow(am, 0.6) * pow(ag, 1.25));
}

auto anf2(double am, double ag) {
    static vectd pars = {4635.9, -21.723, 388.72, -94.873, 1162.3, -4993.5};
    vectd vars = {1., am, ag, am*am, am*ag, ag*ag};
    return dot(pars, vars) / pow(am, 1.55) / pow(ag, 1.095);
}

auto anrh(double am, double ag) {
    return 546.8 / pow(am, 0.6) / pow(ag, 1.12);
}

auto an2(double am, double ag) {
    static vectd pars = {0.18907E-02, -0.50839E-05, 0.16955E-02, -0.19040E-05,
                         0.69806E-03, -0.59716E-01};
    vectd vars = {1., am, ag, am*am, am*ag, ag*ag};
    cout << "an2: " << dot(pars, vars) << endl;
    return dot(pars, vars) * pow(am, 0.5) / pow(ag, 1.02);
}

auto an0(double am, double ag) {
    static vectd pars = {-9.5637, 8.7421, -0.28803, -1.8988, 0.40894, -1.5582,
                        -0.84441E-03, 0.29120E-01, -0.32185, 1.9611};
    for (auto v : pars) cout << v << " ";
    vectd vars = {1., am, ag, am*am, am*ag, ag*ag,
                  am*am*am, am*am*ag, am*ag*ag, ag*ag*ag};
    return dot(pars, vars) * am*am*am / pow(ag, 1.15);
}

// ////////////////////////// //
// ### Complex amplitudes ### //
// ////////////////////////// //
// D2 amplitude (=1)
constexpr auto m_dm2 = 0.6 * 0.06;  // sqrt(get_ddd()/an2(m_ams2, m_gms2));
// D0 amplitude
const auto m_dm0 = 0.6 * exp(m_f0 * ione) *
        sqrt(m_a0 / an0(m_ams0, m_gms0));
// D* amplitude
const auto m_dm1 = 1.4 * exp(m_f1 * ione) * sqrt(m_a1);
// rho(770) amplitude
const auto m_dmrh = 1.1 * exp(m_fr*ione) *
        sqrt(m_ar / anrh(m_amrh, m_gmrh));
// f2 amplitude
const auto m_dmf2 = 1.1 * exp((m_fr+m_ff2)*ione) *
        sqrt(m_af2 / anf2(m_amf2, m_gmf2));
// f0(550) amplitude
const auto m_dmf0 = 1.4 * exp((m_fr+m_ff0)*ione) *
        sqrt(m_af0 / anf0(m_amf0, m_gmf0));
// f(980) amplitude
const auto m_dmbs = exp(m_fb*ione) * sqrt(m_ab / anf0(m_amf098, m_gmf098));
// f(1370) amplitude
const auto m_d3 = 0.15 * exp(m_f3*ione) * sqrt(m_a3 / anf0(m_am3, m_gm3));
// Resonance radius in Blatt-Weiskopf formfactor
constexpr double m_ffr  = 1.6;
constexpr double m_ffdel = 1.0;
// M -> (R -> AB) C. Returns A energy in the R frame
auto eA(double mAsq, double mBsq, double mABsq) {
    return (mAsq + mABsq - mBsq) / (2. * sqrt(mABsq));
}
// M -> (R -> AB) C. Returns A energy in the R frame
auto eB(double mAsq, double mBsq, double mABsq) {
    return (mBsq + mABsq - mAsq) / (2. * sqrt(mABsq));
}
// M -> (R -> AB) C. Returns C energy in the R frame
auto eC(double mMsq, double mCsq, double mABsq) {
    return (mMsq - mABsq - mCsq) / (2. * sqrt(mABsq));
}
// M -> (R -> AB) C. Returns cos of angle between p(A) and p(C) in the R frame
auto MyCosHelAB2(double mMsq, double mAsq, double mBsq, double mCsq,
                 double mABsq, double mBCsq) {
    auto enB = eB(mAsq, mBsq, mABsq);
    auto moBsq = pow(enB, 2) - mBsq;
    if (moBsq < 0) {return make_pair(0., 0.);}
    auto moB = sqrt(moBsq);
    auto enC = eC(mMsq, mCsq, mABsq);
    auto moCsq = pow(enC, 2) - mCsq;
    if (moCsq < 0) {return make_pair(0., 0.);}
    auto moC = sqrt(moCsq);

    return make_pair(
        (mBsq + mCsq + 2.*enB*enC - mBCsq) / (2. * moB * moC), moB * moC
    );
}
// Helicity for a (D pi+) resonance, A = D0, B = pi, C = pi
auto CosHelD0pi2(double mDpipsq, double mpipimsq) {
    return MyCosHelAB2(m_fm2, m_dstm2, m_pm2, m_pm2, mDpipsq, mpipimsq);
}
// A = pi, B = pi, C = D0
auto CosHelpipi2(double mpipisq, double mDpimsq) {
    return MyCosHelAB2(m_fm2, m_pm2, m_pm2, m_dstm2, mpipisq, mDpimsq);
}
// Amplitude for a scalar D0 pi resonance
auto amrs0(double q2) {
    static auto am2 = m_ams0 * m_ams0;
    static auto mom0 = momt(m_ams0, m_pm2, m_dstm2);
    auto q = sqrt(q2);
    auto ppi_frac = momt(q, m_pm2, m_dstm2) / mom0;
    return 1. / compld(q2 - am2, ppi_frac * (m_gms0 * am2 / q));
}
// D*v amplitude
auto amrs1(double q2) {
    // D* veto
//    if (std::fabs(q2 - 4.04) < 0.01) return compld(0, 0);
//    static auto am  = 2.01;
//    static auto gm  = 0.0001;
//    static auto am2 = am * am;
//    static auto ffexpr = [=](double r) {return 1. + r*r;};
//    auto q  = sqrt(q2);
//    // B0 formfactor
//    auto rr_b  = momt(m_fm, q2,  m_pm2) * m_ffr;
//    static auto r0_b  = momt(m_fm, am2, m_pm2) * m_ffr;
//    auto fg = sqrt(ffexpr(r0_b) / ffexpr(rr_b));
//    // resonance formfactor
//    static auto ppi0 = momt(am, m_pm2, m_dstm2);
//    auto ppi  = momt(q,  m_pm2, m_dstm2);
//    auto f = exp((-ppi + ppi0) * m_ffdel);
//    // full Breit-Wigner
//    return f * fg / compld(q2 - am2, (gm * am2 / q) *
//                           pow(ppi / ppi0, 3) * pow(f, 2)) * 3.53438;
    if (std::fabs(q2-4.04) < 0.01) return compld(0, 0);
    const double am    = 2.01;
    const double gm    = 0.0001;
    const double am2   = pow(am, 2);
    const double q     = sqrt(q2);
    // pion energy in the resonance frame calculated with nominal
    // resonance mass
    const double epi0  = (am2+m_pm2-m_dstm2)/(2.*am);
    const double epi   = (q2 +m_pm2-m_dstm2)/(2.*q);
    // pion momentum in the resonance frame calculated with nominal
    // resonance mass
    const double ppi0  = sqrt(pow(epi0, 2)-m_pm2);
    const double ppi   = sqrt(pow(epi , 2)-m_pm2);
    // resonance energy in the B frame calculated with nominal resonance mass
    const double eres0 = (m_fm2+am2-m_pm2)/(2.*m_fm);
    const double eres  = (m_fm2+q2 -m_pm2)/(2.*m_fm);
    // resonance momentum
    const double pres0 = sqrt(pow(eres0, 2)-am2);
    const double pres  = sqrt(pow(eres, 2)-am2);
    // B0 formfactor
    const double rr_b  = pres *m_ffr;
    const double r0_b  = pres0*m_ffr;
    const double fg    = sqrt((1+pow(r0_b, 2))/(1+pow(rr_b, 2)));
    // resonance formfactor
    const double rr_r  = ppi *m_ffdel;
    const double r0_r  = ppi0*m_ffdel;
    const double f     = exp(-rr_r+r0_r);
    // full Breit-Wigner
    const double ar = q2-am2;
    const double gr = (gm*am2/q)*pow(ppi/ppi0, 3)*pow(f, 2);
    return f*fg/(compld(ar, gr)*3.53438);
}
// D*2(2460) amplitude
auto amrs2(double q2, double am, double gm) {
    auto am2 = am * am;
    auto q = sqrt(q2);
    // pion momentum
    auto ppi0 = momt(am, m_pm2, m_dstm2);
    auto ppi  = momt(q,  m_pm2, m_dstm2);
    // resonance momentum
    auto pres0 = momt(m_fm, am2, m_pm2);
    auto pres  = momt(m_fm, q2,  m_pm2);
    // lambda for formfactors
    static auto ffexpr = [=](double r) {return 9. + 2.*r*r + r*r*r*r;};
    // B0 formfactor
    auto fg = sqrt(ffexpr(pres0 * m_ffr) / ffexpr(pres * m_ffr));
    // resonance formfactor
    auto ff = sqrt(ffexpr(ppi0  *  m_ffr) / ffexpr(ppi *  m_ffr));
    // full Breit-Wigner
    return ff * fg /
            compld(q2 - am2, (gm * am2 / q) * pow(ppi / ppi0, 5) * pow(ff, 2));
}
// A scalar pi pi resonance amplitude
auto amrs0h(double q2, double am, double gm) {
    auto am2 = pow(am, 2);
    auto ppi_frac = momt(0.25 * q2, m_pm2) / momt(0.25 * am2, m_pm2);
    return 1./compld(q2 - am2, (gm * am2 / sqrt(q2)) * ppi_frac);
}
// WTF???
auto hwrho(double s) {
    double y = sqrt(1. - 0.07795/s);
    double w = log((1. + y) / (1. - y));
    return make_pair(
                w * y * M_PI_2,
                (0.5 * w * (1. - y*y) / y + 1.) / s * M_PI_2
                );
}
// rho(770) and rho() amplitude
auto amrs1h(double q2, double am, double gm) {
    auto am2  = pow(am, 2);
    // pion momentum
    auto ppi0 = momt(0.25 * am2, m_pm2);
    auto ppi  = momt(0.25 * q2,  m_pm2);
    // resonance energy
    auto ar  = q2 - am2;
    auto gr  = (gm * am2 / sqrt(q2)) * pow(ppi / ppi0, 3);
    auto hw  = hwrho(q2).first;
    auto tmp = hwrho(am2);
    auto hw0 = tmp.first;
    auto hwd0 = tmp.second;
    auto dm  = ((hw - hw0) * pow(ppi, 2) - ar * hwd0 * pow(ppi0, 2))
               * am2 * gm / pow(ppi0, 3);
    return 1. / compld(ar - dm, gr);
}
// f2 tensor pi pi resonance amplitude
auto amrs2h(double q2) {
    static auto am2 = pow(m_amf2, 2);
    // pion momentum
    auto ppi_frac = momt(0.25 * q2,  m_pm2) / momt(0.25 * am2, m_pm2);
    return 1. / compld(q2 - am2, (m_gmf2 * am2 / sqrt(q2)) * pow(ppi_frac, 5));
}
// Rho-omega interference amplitude
auto amrs1hm(double q2) {
    static double am2 = pow(0.78257, 2);
    return 1. / compld(q2 - am2, sqrt(q2) * 0.00849);
}
// Angular factor for vectors
inline auto ang1(double h, double pq) {return pq * h;}
// Angular factor for tensors
inline auto ang2(double h, double pq) {return pq*pq * (h*h - 1./3.);}
// D* amplitude
auto DstarAmp(double qp, double h2, double pq2) {
    return ang1(h2, pq2) * amrs1(qp);
}
// D** amplitude
auto DdstarAmp(double qp, double h2, double pq2) {
    return ang2(h2, pq2)*amrs2(qp, m_ams2, m_gms2);
}
// D0* amplitude
auto DzeroAmp(double qp) {return amrs0(qp);}
// rho(770) amplitude
auto rhoAmp(double q3, double h3, double pq3) {
    return ang1(h3, pq3) * amrs1h(q3, m_amrh, m_gmrh);
}
// omega amplitude
auto omegaAmp(double q3, double h3, double pq3) {
    return ang1(h3, pq3) * amrs1hm(q3) * m_amom * exp(2.93 * ione);
}
// rho' amplitude
auto rhopAmp(double q3, double h3, double pq3) {
    return ang1(h3, pq3) * amrs1h(q3, m_amrh1, m_gmrh1) * m_rhop_amp;
}
// rho(770) + omega + rho' amplitude
auto FullRhoAmp(double q3, double h3, double pq3) {
    return ang1(h3, pq3) * (amrs1h(q3, m_amrh, m_gmrh)
                         + amrs1hm(q3) * m_amom * exp(1.99 * ione)
                         + amrs1h(q3, m_amrh1, m_gmrh1) * m_rhop_amp);
}
// f2 amplitude
auto f2Amp(double q3, double h3, double pq3) {
    return ang2(h3, pq3) * amrs2h(q3);
}
// f0 amplitude
auto f0Amp(double q3) {
    return amrs0h(q3, m_amf0, m_gmf0);
}
// f(980) amplitude
auto f980Amp(double q3) {
    return amrs0h(q3, m_amf098, m_gmf098);
}
// f(1370) amplitude
std::complex<double> f1370Amp(double q3) {
    return amrs0h(q3, m_am3, m_gm3);
}

KAmp::KuzBtoDpipiAmp() :
    AbsDalitzModel(m_fm, m_dstm, m_pm, m_pm),
    AbsSymDalitzModel(m_fm, m_dstm, m_pm) {
    SetModelTitle("Belle amplitude (Kuzmin)");
    SetResNames({"D*", "D2", "D0", "rho", "omega", "rho'", "f2", "f0(550)",
                 "f(980)", "f(1370)"});
    SetAmpNames({"D*", "D2", "D0", "rho", "f2", "f0(550)",
                 "f(980)", "f(1370)"});
    SetABaxis("m^{2}(D\\pi^{+})\\ (GeV^{2}/c^{4})");
    SetACaxis("m^{2}(D\\pi^{-})\\ (GeV^{2}/c^{4})");
    SetBCaxis("m^{2}(\\pi^{+}\\pi^{-})\\ (GeV^{2}/c^{4})");
    init();
    PrintSummary();
}

compld KAmp::GetResVal(double mABsq, double mACsq, const int resnum) const {
    double pq  = 0;
    double hel = 0;
    const double mBCsq = m3sq(mABsq, mACsq);
    if (resnum < 2) {
        auto tmp = CosHelD0pi2(mABsq, mBCsq);  // D pi resonance
        hel = tmp.first, pq = tmp.second;
    } else if (resnum < 5) {
        auto tmp = CosHelpipi2(mBCsq, mABsq);  // pi pi resonance
        hel = tmp.first, pq = tmp.second;
    }
    switch (resnum) {
    case 0: return DstarAmp(mABsq, hel, pq);  // D*
    case 1: return DdstarAmp(mABsq, hel, pq);  // D2
    case 2: return DzeroAmp(mABsq);  // D0
    case 3: return FullRhoAmp(mBCsq, hel, pq);  // rho + omega + rho'
    case 4: return f2Amp(mBCsq, hel, pq);  // f2
    case 5: return f0Amp(mBCsq);  // f0
    case 6: return f980Amp(mBCsq);  // f980
    case 7: return f1370Amp(mBCsq);  // f1370
    default: cerr << "KAmp: wrong resnum " << resnum << endl;
    }
    return 0;
}

void KAmp::GetResVals(vectcd* resv, double mABsq, double mACsq) const {
    auto mBCsq = m3sq(mABsq, mACsq);
    double h2, pq2;
    double h3, pq3;
    auto tmp2 = CosHelD0pi2(mABsq, mBCsq);
    h2 = tmp2.first, pq2 = tmp2.second;
    auto tmp3 = CosHelpipi2(mBCsq, mACsq);
    h3 = tmp3.first, pq3 = tmp3.second;
    *resv = {
        DstarAmp(mABsq, h2, pq2),  // D* amplitude
        DdstarAmp(mABsq, h2, pq2),  // D** amplitude
        DzeroAmp(mABsq),  // D0* amplitude
        rhoAmp(mBCsq, h3, pq3),  // rho(770) amplitude
        omegaAmp(mBCsq, h3, pq3),  // omega amplitude
        rhopAmp(mBCsq, h3, pq3),  // rho' amplitude
        f2Amp(mBCsq, h3, pq3),  // f2 amplitude
        f0Amp(mBCsq),  // f0 amplitude
        f980Amp(mBCsq),  // f(980) amplitude
        f1370Amp(mBCsq)  // f(1370) amplitude
    };
}

void KAmp::init(void) {
    SetCoefficients({m_dm1, m_dm2, m_dm0, m_dmrh, m_dmf2,
                     m_dmf0, m_dmbs, m_d3});
    SetAmpSignature({1, 1, 1, 3, 1, 1, 1, 1});
}

double KAmp::Norm(void) const {return NormWithCache();}

void KAmp::SetParams(const vectd& pars) {
    // parameters -> amplitudes and phases
    if (pars.size() != 2 * AmpNum()) {
        cerr << "Wrong pars size: " << pars.size() << " (" << 2*AmpNum()
             << " expected)" << endl;
        return;
    }
    vectcd coeffs;
    for (auto i = 0u; i < AmpNum(); i++) {
        auto amp = pars[2*i];
        auto pha = pars[2*i+1];
        coeffs.emplace_back(amp * exp(ione * pha));
    }
    SetCoefficients(coeffs);
}

auto cmplStr(const compld& x) {
    return to_string(abs(x)) + " * e^(" + to_string(arg(x)) + ")";
}

void KAmp::PrintSummary(void) const {
    cout << "*** The Kuzmin amplitude ***" << endl;
    cout << "Masses: B0 " << mM() << ", D0 " << mA()
         << ", pi+ " << mB() << endl;
    const auto& coefs = GetCoefficients();
    cout << "D*      amp: " << cmplStr(coefs[0])  << endl;
    cout << "D2      amp: " << cmplStr(coefs[1])  << endl;
    cout << "D0      amp: " << cmplStr(coefs[2])  << endl;
    cout << "rho     amp: " << cmplStr(coefs[3]) << ", anrh: "
         << anrh(m_amrh, m_gmrh) << endl;
    cout << "omega   amp: " << cmplStr(coefs[3]) << endl;
    cout << "rho'    amp: " << cmplStr(coefs[3]) << endl;
    cout << "f2      amp: " << cmplStr(coefs[4]) << ", anf2: "
         << anf2(m_amf2, m_gmf2) << endl;
    cout << "f0(550) amp: " << cmplStr(coefs[6]) << ", anf0: "
         << anf0(m_amf0, m_gmf0) << endl;
    cout << "f(980)  amp: " << cmplStr(coefs[5]) << endl;
    cout << "f(1370) amp: " << cmplStr(coefs[7])   << ", anf0: "
         << anf0(m_am3, m_gm3)  << endl;
}
