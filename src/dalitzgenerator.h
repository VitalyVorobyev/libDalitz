#ifndef DALITZGENERATOR_H
#define DALITZGENERATOR_H

#include <vector>
#include "dalitzmodel.h"
#include "randomdalitzpoint.h"

/// \brief Class for generation of Dalitz plot distribution.
/// Can generate events for any amplitude defined with DalitzModel class

class DalitzGenerator : public RandomDalitzPoint{
public:
  DalitzGenerator(const DalitzModel& _dm);
  int Generate(const int NEv,std::vector<double>& mABv,std::vector<double>& mACv);

  void SetMaxTries(const long& p) {m_ntries = p; return;}
  long GetMaxTries(void) const {return m_ntries;}
  void SetNMajCounts(const int p){m_maj_counts = p; return;}
  void SetMajorant(const double& p){m_maj = p; return;}

private:
  double CalcMajorant(void);
  long m_ntries;
  int m_maj_counts;
  double m_maj;
  DalitzModel* dm;
};

#endif // DALITZGENERATOR_H
