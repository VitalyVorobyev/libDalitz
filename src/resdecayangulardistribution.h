#ifndef RESDECAYANGULARDISTRIBUTION_H
#define RESDECAYANGULARDISTRIBUTION_H

/// Class for computation of angular pdf for intermidiate resonance decay in a three-body decay.
/// Final state particles are assumed to be scalars.
/// Resonance may be scalar, vector or tensor.

class ResDecayAngularDistribution{
public:
  ResDecayAngularDistribution(const int spin, const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres);
  double operator()(const double& mACsq,const double& mBCsq) const;

  void Set_mR(const double& mr) {m_mre_sq = mr*mr; return;}

  double mMotSq(void) const {return m_mmo_sq;}
  double mChASq(void) const {return m_mca_sq;}
  double mChBSq(void) const {return m_mcb_sq;}
  double mChCSq(void) const {return m_mcc_sq;}
  double mResSq(void) const {return m_mre_sq;}

private:
  int m_spin;
  double m_mmo_sq;
  double m_mca_sq;
  double m_mcb_sq;
  double m_mcc_sq;
  double m_mre_sq;

  double m_pc_var1;
  double m_pc_var2;
  double m_pc_var3;
  double m_pc_var4;
  double m_pc_var5;

  void make_precalc_one(void);
  void make_precalc_all(void);
};

#endif // RESDECAYANGULARDISTRIBUTION_H
