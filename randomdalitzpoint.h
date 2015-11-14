#ifndef RANDOMDALITZPOINT_H
#define RANDOMDALITZPOINT_H

#include "dalitzphasespace.h"
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>

class RandomDalitzPoint : public DalitzPhaseSpace{
public:
  RandomDalitzPoint(const DalitzPhaseSpace& phsp);
  RandomDalitzPoint(const double& mmo, const double& mca, const double& mcb, const double& mcc);
  ~RandomDalitzPoint();

  void GetPoint(double& mAB, double& mAC);
  void SetSeed(const int seed);
  unsigned GetSeed(void) const {return m_seed;}

private:
  void init();
  std::uniform_real_distribution<double>* unifAB;
  std::uniform_real_distribution<double>* unifAC;
  std::default_random_engine re;
  unsigned m_seed;
};

#endif // RANDOMDALITZPOINT_H
