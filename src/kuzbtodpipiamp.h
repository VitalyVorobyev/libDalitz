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

/**
 * \brief The KuzBtoDpipiAmp class. This class reproducese Alex Kuzmin's
 * Fortran code describing the B0 -> D0 pi+ pi- decay amplitude model
 * obtained with Belle data for Ref!!!
 **/
class KuzBtoDpipiAmp : public AbsSymDalitzModel {
 public:
    KuzBtoDpipiAmp();
    /// Vector of Res_amp (rho, omega and rho' combined)
    void GetResVals(std::vector<std::complex<double>>* resv,
                    const double& mABsq, const double& mACsq) const;
    ///
    std::complex<double> GetResVal(const double& mABsq, const double& mACsq,
                                   const int resnum) const;
    /// Set parameters
    void SetParams(const std::vector<double>& pars);
    /// Get normalization
    double Norm(void) const;

 private:
    void PrintSummary(void);
    void init(void);
    /// Cos of helisity angle for resonance R -> AB
    static double CosHelAB(const double& mMo, const double& mA,
                           const double& mB, const double& mC,
                           const double& mABsq);
    /// Helicity for a (D pi+) resonance
    static double CosHelD0pi(const double& mDpipsq, const double& mDpimsq,
                             double* pq);
    static double CosHelD0pi2(const double& mDpipsq, const double& mpipimsq,
                              double* pq);
    /// Helicity for a (pi+ pi-) resonance
    static double CosHelpipi(const double& mpipisq, const double& mDpipsq,
                             double* pq);
    static double CosHelpipi2(const double& mpipisq, const double& mDpimsq,
                              double* pq);
    /// M -> (R -> AB) C. Returns A energy in the R frame
    static double eA(const double& mAsq, const double& mBsq,
                     const double& mABsq);
    /// M -> (R -> AB) C. Returns B energy in the R frame
    static double eB(const double& mAsq, const double& mBsq,
                     const double& mABsq);
    /// M -> (R -> AB) C. Returns C energy in the R frame
    static double eC(const double& mMsq, const double& mCsq,
                     const double& mABsq);
    /// M -> (R -> AB) C. Returns C energy in the M frame
    static double eC0(const double& mMsq, const double& mCsq,
                      const double& mABsq);
    /// M -> (R -> AB) C. Returns cos of angle between p(A) and p(C)
    /// in the R frame
    static double MyCosHelAB(const double& mMsq, const double& mAsq,
                             const double& mBsq, const double& mCsq,
                             const double& mABsq, const double& mACsq,
                             double* pq);
    static double MyCosHelAB2(const double& mMsq, const double& mAsq,
                              const double& mBsq, const double& mCsq,
                              const double& mABsq, const double& mBCsq,
                              double* pq);
    /// M -> (R -> AB) C. Returns Dot of p(A) and p(C) 3-vectors in the R frame
    static double DotAC(const double& mMsq, const double& mAsq,
                        const double& mBsq, const double& mCsq,
                        const double& mABsq, const double& mACsq);
    /// Summ of all amplitudes (???)
    static double get_ddd(void);
    /// Complex amplitude for D0
    static std::complex<double> get_dm0(void);
    /// Complex amplitude for D*
    static std::complex<double> get_dm1(void);
    /// Complex amplitude for D2
    static double get_dm2(void);
    /// Complex amplitude for f(980)
    static std::complex<double> get_dmbs(void);
    /// Complex amplitude for f(1370)
    static std::complex<double> get_d3(void);
    /// Complex amplitude for rho(770)
    static std::complex<double> get_dmrh(void);
    /// Complex amplitude for f(550)
    static std::complex<double> get_dmf0(void);
    /// Complex amplitude for f2
    static std::complex<double> get_dmf2(void);

    /// D* amplitude
    static std::complex<double> DstarAmp(const double& qp, const double& h2,
                                         const double& pq2);
    /// D** amplitude
    static std::complex<double> DdstarAmp(const double& qp, const double& h2,
                                          const double& pq2);
    /// D0* amplitude
    static std::complex<double> DzeroAmp(const double& qp);
    /// rho(770) amplitude
    static std::complex<double> rhoAmp(const double& q3, const double& h3,
                                       const double& pq3);
    /// omega amplitude
    static std::complex<double> omegaAmp(const double& q3, const double& h3,
                                         const double& pq3);
    /// rho' amplitude
    static std::complex<double> rhopAmp(const double& q3, const double& h3,
                                        const double& pq3);
    /// rho(770) + omega + rho' amplitude
    static std::complex<double> FullRhoAmp(const double& q3, const double& h3,
                                           const double& pq3);
    /// f2 amplitude
    static std::complex<double> f2Amp(const double& q3, const double& h3,
                                      const double& pq3);
    /// f0 amplitude
    static std::complex<double> f0Amp(const double& q3);
    /// f(980) amplitude
    static std::complex<double> f980Amp(const double& q3);
    /// f(1370) amplitude
    static std::complex<double> f1370Amp(const double& q3);

    /// Amplitude for a scalar D0 pi resonance
    static std::complex<double> amrs0(const double& q2, const double& am,
                                      const double& gm);
    /// Amplitude for a vector D0 pi resonance (only D* case)
    static std::complex<double> amrs1(const double& q2);
    /// Amplitude for a tensor D0 pi resonance
    static std::complex<double> amrs2(const double& q2, const double& am,
                                      const double& gm);
    /// Amplitude for a scalar pi pi resonance
    static std::complex<double> amrs0h(const double& q2, const double& am,
                                       const double& gm);
    /// Amplitude for a vector pi pi resonance
    static std::complex<double> amrs1h(const double& q2, const double& am,
                                       const double& gm);
    /// Amplitude for a tensor pi pi resonance
    static std::complex<double> amrs2h(const double& q2, const double& am,
                                       const double& gm);
    /// Amplitude with rho-omega interference
    static std::complex<double> amrs1hm(const double& q2);
    ///
    static double hwrho(const double& s, double* hwd1);

    static double anf0(const double& am, const double& ag);
    static double anf2(const double& am, const double& ag);
    static double anrh(const double& am, const double& ag);
    static double an2(const double& am, const double& ag);
    static double an0(const double& am, const double& ag);

    static double ang1(const double& h, const double& pq);
    static double ang2(const double& h, const double& pq);

    /// par 1: 'raw' amplitude for D0
    static const double m_a0;
    /// par 3: phase for D0
    static const double m_f0;
    /// par 2: 'raw' amplitude for D*
    static const double m_a1;
    /// par 4: phase for D*
    static const double m_f1;
    /// par 10: 'raw' amplitude for f(980)
    static const double m_ab;
    /// par 11: phase for f(980)
    static const double m_fb;
    /// par 24: 'raw' amplitude for f(1370)
    static const double m_a3;
    /// par 25: phase for f(1370)
    static const double m_f3;
    /// par 14: 'raw' amplitude for rho(770)
    static const double m_ar;
    /// par 15: phase for rho(770)
    static const double m_fr;
    /// par 18: 'raw' amplitude for f0(550)
    static const double m_af0;
    /// par 19: phase for f0(550) relative to the rho(770) phase
    static const double m_ff0;
    /// par 22: 'raw' amplitude for f2
    static const double m_af2;
    /// par 23: phase for f2 relative to the rho(770) phase
    static const double m_ff2;

    /// B0 mass
    static const double m_fm;
    /// B0 mass sq
    static const double m_fm2;
    /// pi+ mass
    static const double m_pm;
    /// pi+ mass sq
    static const double m_pm2;
    /// D0 mass
    static const double m_dstm;
    /// D0 mass sq
    static const double m_dstm2;
    /// par 16: f0 mass
    static const double m_amf0;
    /// par 17: f0 width
    static const double m_gmf0;
    /// par 20: f2 mass
    static const double m_amf2;
    /// par 21: f2 width
    static const double m_gmf2;
    /// f(980) mass
    static const double m_amf098;
    /// f(980) width
    static const double m_gmf098;
    /// f(1370) mass
    static const double m_am3;
    /// f(1370) width
    static const double m_gm3;
    /// D*0 mass
    static const double m_ams0;
    /// D*0 width
    static const double m_gms0;
    /// D2 mass
    static const double m_ams2;
    /// D2 width
    static const double m_gms2;
    /// par 12: rho(770) mass
    static const double m_amrh;
    /// par 13: rho(770) width
    static const double m_gmrh;
    /// Mass rho'
    static const double m_amrh1;
    /// Width rho'
    static const double m_gmrh1;

    /// Resonance radius in Blatt-Weiskopf formfactor
    static const double m_ffr;
    ///
    static const double m_ffdel;

    /// par 30: omega amplitude
    static const double m_amom;
    /// amplitudes sum
    static const double m_ddd;

    /// D2 amplitude (=1)
    static const std::complex<double> m_dm2;
    /// D0 amplitude
    static const std::complex<double> m_dm0;
    /// D* amplitude
    static const std::complex<double> m_dm1;
    /// rho(770) amplitude
    static const std::complex<double> m_dmrh;
    /// f2 amplitude
    static const std::complex<double> m_dmf2;
    /// f0(550) amplitude
    static const std::complex<double> m_dmf0;
    /// f(980) amplitude
    static const std::complex<double> m_dmbs;
    /// f(1370) amplitude
    static const std::complex<double> m_d3;

    static const double m_fm2_2pm2_md2;
};

#endif  // SRC_KUZBTODPIPIAMP_H_
