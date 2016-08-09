#ifndef RANDOMDALITZPOINT_H
#define RANDOMDALITZPOINT_H

#include "dalitzphasespace.h"
#include "dstrip.h"
//#include <algorithm>    // std::shuffle
#include <random> // std::default_random_engine

typedef std::uniform_real_distribution<double> unif;
typedef std::normal_distribution<double>       gaus;
typedef unsigned long long                     U64;
typedef const double                           cdouble;

///
/// \brief The RandomDalitzPoint class
///
class RandomDalitzPoint : public DalitzPhaseSpace{
public:
  RandomDalitzPoint(const DalitzPhaseSpace* phsp);
  RandomDalitzPoint(cdouble& mmo, cdouble& mca, cdouble& mcb, cdouble& mcc);
  ~RandomDalitzPoint();

  int GetPoint(     double& mABsq, double& mACsq) const;
  int GetUnconstrainedPoint(double& mABsq, double& mACsq) const;
  int GetStripPoint(double& mABsq, double& mACsq, const DStrip* shape) const;
  int GetGaussPoint(double& mABsq, double& mACsq, double& mBCsq, const DStrip* shape) const;
  int GetStripPoint(double& mABsq, double& mACsq, const DStrip* shape1, const DStrip* shape2) const;
  int GetGaussPoint(double& mABsq, double& mACsq, double& mBCsq, const DStrip* shape1, const DStrip* shape2) const;

  int GetPoints(     vectd& mABsqV,vectd& mACsqV,const int N) const;
  int GetStripPoints(vectd& mABsqV,vectd& mACsqV,const int N, const DStrip* shape) const;
  int GetGaussPoints(vectd& mABsqV,vectd& mACsqV,vectd& mBCsqV, const int N, const DStrip* shape) const;
  int GetStripPoints(vectd& mABsqV,vectd& mACsqV,const int N, const DStrip* shape1, const DStrip* shape2) const;
  int GetGaussPoints(vectd& mABsqV,vectd& mACsqV,vectd& mBCsqV,const int N, const DStrip* shape1, const DStrip* shape2) const;

  // * Static methods * //
  static void SetSeed(const unsigned seed);
  static unsigned GetSeed(void) {return m_seed;}

private:
  void init();

  int GetUnifPoint(double& mABsq,double& mACsq,const DStrip* shape, unif& dist) const;
  int GetGausPoint1(double& mABsq,double& mACsq, double& mBCsq,const DStrip* shape, gaus& dist) const;

  int GetUnifPoint(double& mABsq,double& mACsq,const DStrip* shape1, unif& dist1,const DStrip* shape2, unif& dist2) const;
  int GetGausPoint2(double& mABsq,double& mACsq, double& mBCsq,const DStrip* shape1, gaus& dist1,const DStrip* shape2, gaus& dist2) const;

  U64 m_max_tries;
  cdouble mABsqMin;
  cdouble mABsqMax;
  cdouble mACsqMin;
  cdouble mACsqMax;

//  unif* unifAB;
//  unif* unifAC;
  static std::default_random_engine re;
  static unsigned m_seed;


};

#endif // RANDOMDALITZPOINT_H
