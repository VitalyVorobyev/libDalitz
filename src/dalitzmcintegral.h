#ifndef DALITZMCINTEGRAL_H
#define DALITZMCINTEGRAL_H

#include "dalitzmodel.h"
#include "randomdalitzpoint.h"

class DalitzMCIntegral : public RandomDalitzPoint{
public:
  DalitzMCIntegral(const DalitzModel& _dm);
  double GetIntegral(const long& nc = 0);

  void SetNCounts(const long& p){m_ncounts = p; return;}
  long GetNCounts(void) const {return m_ncounts;}
private:
  DalitzModel* dm;
  long m_ncounts;
  double m_int;
  double m_int_err;
};

#endif // DALITZMCINTEGRAL_H
