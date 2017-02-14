/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzstrip.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DALITZSTRIP_H_
#define SRC_DALITZSTRIP_H_

#include "./dalitzveto.h"

class DalitzStrip : public DalitzVeto {
 public:
    DalitzStrip();

 private:
    int m_type;
};

#endif  // SRC_DALITZSTRIP_H_
