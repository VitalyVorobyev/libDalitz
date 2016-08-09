#ifndef KUZBTODPIPIAMP_H
#define KUZBTODPIPIAMP_H

#include <vector>
#include <complex>
#include <string>

#include "absdalitzmodel.h"

///
/// \brief The KuzBtoDpipiAmp class. This class reproducese Alex Kuzmin's Fortran code
/// describing the B0 -> D0 pi+ pi- decay amplitude model obtained with Belle data
/// for Ref!!!
///
class KuzBtoDpipiAmp : public AbsDalitzModel{
public:
  KuzBtoDpipiAmp();
  /// Vector of Res_amp (rho, omega and rho' combined)
  void GetResVals(vectcd& resv,cdouble& mABsq,cdouble& mACsq) const;
  ///
  compld GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const;
  /// Set parameters
  void SetParams(const vectd& pars);
  /// Get normalization
  double Norm(void) const;
private:
  void PrintSummary(void);
  void init(void);

  static double GetmBCsq(cdouble& mABsq,cdouble& mACsq);
  /// Cos of helisity angle for resonance R -> AB
  static double CosHelAB(cdouble& mMo, cdouble& mA, cdouble& mB, cdouble& mC, cdouble& mABsq);
  /// Helicity for a (D pi+) resonance
  static double CosHelD0pi( cdouble& mDpipsq, cdouble& mDpimsq,  double& pq);
  static double CosHelD0pi2(cdouble& mDpipsq, cdouble& mpipimsq, double& pq);
  /// Helicity for a (pi+ pi-) resonance
  static double CosHelpipi( cdouble& mpipisq, cdouble& mDpipsq, double& pq);
  static double CosHelpipi2(cdouble& mpipisq, cdouble& mDpimsq, double& pq);

  /// M -> (R -> AB) C. Returns A energy in the R frame
  static double eA(cdouble& mAsq, cdouble& mBsq, cdouble& mABsq);
  /// M -> (R -> AB) C. Returns B energy in the R frame
  static double eB(cdouble& mAsq, cdouble& mBsq, cdouble& mABsq);
  /// M -> (R -> AB) C. Returns C energy in the R frame
  static double eC(cdouble& mMsq, cdouble& mCsq, cdouble& mABsq);
  /// M -> (R -> AB) C. Returns C energy in the M frame
  static double eC0(cdouble& mMsq, cdouble& mCsq, cdouble& mABsq);
  /// M -> (R -> AB) C. Returns cos of angle between p(A) and p(C) in the R frame
  static double MyCosHelAB(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq, double& pq);
  static double MyCosHelAB2(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mBCsq, double& pq);
  /// M -> (R -> AB) C. Returns Dot of p(A) and p(C) 3-vectors in the R frame
  static double DotAC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq);

  /// Summ of all amplitudes (???)
  static double get_ddd(void);

  /// Complex amplitude for D0
  static compld get_dm0(void);
  /// Complex amplitude for D*
  static compld get_dm1(void);
  /// Complex amplitude for D2
  static double get_dm2(void);
  /// Complex amplitude for f(980)
  static compld get_dmbs(void);
  /// Complex amplitude for f(1370)
  static compld get_d3(void);
  /// Complex amplitude for rho(770)
  static compld get_dmrh(void);
  /// Complex amplitude for f(550)
  static compld get_dmf0(void);
  /// Complex amplitude for f2
  static compld get_dmf2(void);

//  static double mom_mult(cdouble& q2);

//  compld btodpp(cdouble& qp, cdouble& qn) const;

  /// D* amplitude
  static compld DstarAmp(cdouble& qp,cdouble& h2,cdouble& pq2);
  /// D** amplitude
  static compld DdstarAmp(cdouble& qp,cdouble& h2,cdouble& pq2);
  /// D0* amplitude
  static compld DzeroAmp(cdouble& qp);
  /// rho(770) amplitude
  static compld rhoAmp(cdouble& q3,cdouble& h3,cdouble& pq3);
  /// omega amplitude
  static compld omegaAmp(cdouble& q3,cdouble& h3,cdouble& pq3);
  /// rho' amplitude
  static compld rhopAmp(cdouble& q3,cdouble& h3,cdouble& pq3);
  /// rho(770) + omega + rho' amplitude
  static compld FullRhoAmp(cdouble& q3,cdouble& h3,cdouble& pq3);
  /// f2 amplitude
  static compld f2Amp(cdouble& q3,cdouble& h3,cdouble& pq3);
  /// f0 amplitude
  static compld f0Amp(cdouble& q3);
  /// f(980) amplitude
  static compld f980Amp(cdouble& q3);
  /// f(1370) amplitude
  static compld f1370Amp(cdouble& q3);

  /// Amplitude for a scalar D0 pi resonance
  static compld amrs0(cdouble& q2,cdouble& am,cdouble& gm);
  /// Amplitude for a vector D0 pi resonance (only D* case)
  static compld amrs1(cdouble& q2);
  /// Amplitude for a tensor D0 pi resonance
  static compld amrs2(cdouble& q2,cdouble& am,cdouble& gm);
  /// Amplitude for a scalar pi pi resonance
  static compld amrs0h(cdouble& q2,cdouble& am,cdouble& gm);
  /// Amplitude for a vector pi pi resonance
  static compld amrs1h(cdouble& q2,cdouble& am,cdouble& gm);
  /// Amplitude for a tensor pi pi resonance
  static compld amrs2h(cdouble& q2,cdouble& am,cdouble& gm);
  /// Amplitude with rho-omega interference
  static compld amrs1hm(cdouble& q2);
  ///
  static double hwrho(cdouble& s, double& hwd1);

  static double anf0(cdouble& am, cdouble& ag);
  static double anf2(cdouble& am, cdouble& ag);
  static double anrh(cdouble& am, cdouble& ag);
  static double an2( cdouble& am, cdouble& ag);
  static double an0( cdouble& am, cdouble& ag);

  static double ang1(cdouble& h,cdouble& pq);
  static double ang2(cdouble& h,cdouble& pq);

  /// par 1: 'raw' amplitude for D0
  static cdouble m_a0;
  /// par 3: phase for D0
  static cdouble m_f0;
  /// par 2: 'raw' amplitude for D*
  static cdouble m_a1;
  /// par 4: phase for D*
  static cdouble m_f1;
  /// par 10: 'raw' amplitude for f(980)
  static cdouble m_ab;
  /// par 11: phase for f(980)
  static cdouble m_fb;
  /// par 24: 'raw' amplitude for f(1370)
  static cdouble m_a3;
  /// par 25: phase for f(1370)
  static cdouble m_f3;
  /// par 14: 'raw' amplitude for rho(770)
  static cdouble m_ar;
  /// par 15: phase for rho(770)
  static cdouble m_fr;
  /// par 18: 'raw' amplitude for f0(550)
  static cdouble m_af0;
  /// par 19: phase for f0(550) relative to the rho(770) phase
  static cdouble m_ff0;
  /// par 22: 'raw' amplitude for f2
  static cdouble m_af2;
  /// par 23: phase for f2 relative to the rho(770) phase
  static cdouble m_ff2;

  /// B0 mass
  static cdouble m_fm;
  /// B0 mass sq
  static cdouble m_fm2;
  /// pi+ mass
  static cdouble m_pm;
  /// pi+ mass sq
  static cdouble m_pm2;
  /// D0 mass
  static cdouble m_dstm;
  /// D0 mass sq
  static cdouble m_dstm2;
  /// par 16: f0 mass
  static cdouble m_amf0;
  /// par 17: f0 width
  static cdouble m_gmf0;
  /// par 20: f2 mass
  static cdouble m_amf2;
  /// par 21: f2 width
  static cdouble m_gmf2;
  /// f(980) mass
  static cdouble m_amf098;
  /// f(980) width
  static cdouble m_gmf098;
  /// f(1370) mass
  static cdouble m_am3;
  /// f(1370) width
  static cdouble m_gm3;
  /// D*0 mass
  static cdouble m_ams0;
  /// D*0 width
  static cdouble m_gms0;
  /// D2 mass
  static cdouble m_ams2;
  /// D2 width
  static cdouble m_gms2;
  /// par 12: rho(770) mass
  static cdouble m_amrh;
  /// par 13: rho(770) width
  static cdouble m_gmrh;
  /// Mass rho'
  static cdouble m_amrh1;
  /// Width rho'
  static cdouble m_gmrh1;

  /// Resonance radius in Blatt-Weiskopf formfactor
  static cdouble m_ffr;
  ///
  static cdouble m_ffdel;

  /// par 30: omega amplitude
  static cdouble m_amom;
  /// amplitudes sum
  static cdouble m_ddd;

  /// D2 amplitude (=1)
  static const compld m_dm2;
  /// D0 amplitude
  static const compld m_dm0;
  /// D* amplitude
  static const compld m_dm1;
  /// rho(770) amplitude
  static const compld m_dmrh;
  /// f2 amplitude
  static const compld m_dmf2;
  /// f0(550) amplitude
  static const compld m_dmf0;
  /// f(980) amplitude
  static const compld m_dmbs;
  /// f(1370) amplitude
  static const compld m_d3;

  static cdouble m_fm2_2pm2_md2;
};

#endif // KUZBTODPIPIAMP_H
