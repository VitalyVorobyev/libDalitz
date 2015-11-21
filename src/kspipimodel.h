#ifndef KspipiModel_H
#define KspipiModel_H

#include "symdalitzmodel.h"

#include <math.h>
#include <vector>

/// \brief Implementation of D0 -> Ks0 pi+ pi- Belle 2010 decay model.
/// Parameters are taken from A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010.
/// See also S. Kopp et al. (CLEO Collaboration) Phys. Rev. D 63, 092001 – Published 9 April 2001

class KspipiModel : public SymDalitzModel{
public:
  KspipiModel(void);
  KspipiModel(const double& md, const double& mks, const double& mpi);
private:
};

#endif // KspipiModel_H
