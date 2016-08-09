#ifndef SMALLDPBIN_H
#define SMALLDPBIN_H

#include "dalitzmodel.h"
#include "consts.h"

///
/// \brief The SmallDPBin class
///
class SmallDPBin{
public:
  SmallDPBin(const double& mAB = 0, const double& mAC = 0);
  SmallDPBin(DalitzModel* model, const double& mAB = 0, const double& mAC = 0);
  SmallDPBin(const SmallDPBin& x);
  SmallDPBin& operator=(const SmallDPBin& x) = default;

  int SetModel(DalitzModel* model);
  int SetGridSize(const int gsize);
  int SetGridStep(const double& gstep);
  int SetPoint(const double& mAB, const double& mAC);
  int SetBin(const int i, const int j);
  int SetNorm(const double& norm);

  bool IsInBin(const double& mAB, const double& mAC) const;
  double Phase(void)  const;
  double Amp(void)    const;
  double Ampb(void)   const;
  double P(void)      const;
  double Pb(void)     const;
  compld CAmp(void)   const;
  compld CAmpb(void)  const;
  int    GSize(void)  const;
  double Weight(void) const;
  int GetCurrentPoint(double& mAB, double& mAC) const;

  int GetBin(const double& x) const;
  double GetVal(const int bin) const;

private:
  int Calc(void);

  DalitzModel* m_model;

  /// Amplitude value for m_mAB and m_mAC
  compld m_val;
  /// Amplitude value for m_mAC and m_mAB
  compld m_valb;
  /// Statistical weight of the bin
  double m_wght;
  /// Phase difference
  double m_dph;
  /// Norm amp
  double m_p;
  /// Norm anti-amp
  double m_pb;

  /// Size of the grid cell
  double m_h;
  /// Number of the grid cells
  int    m_gsize;
  /// Normalization coeffitient for amplitude
  double m_norm;
  /// Range of the grid
  double m_dpsize;
  /// Min value of the Dalitz variable
  double m_dpmin;

  /// Current var x
  double m_mAB;
  /// Current var y
  double m_mAC;
  /// Current bin x
  int m_ABbin;
  /// Current bin y
  int m_ACbin;
};

#endif // SMALLDPBIN_H
