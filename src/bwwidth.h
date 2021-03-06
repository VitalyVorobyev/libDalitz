/** Copyright 2017 Vitaly Vorobyev
 ** @file bwwidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BWWIDTH_H_
#define SRC_BWWIDTH_H_

#include "./absvarwidth.h"
#include "./blattweisskopf.h"

/**
 * @brief The BWWidth class. Eq.(10) in Phys. Rev D92, 032002 (2015)
 */
class BWWidth : public AbsVarWidth {
 public:
    BWWidth(const double& G0, const double& m,
            const double& p0, const int mom);

    double operator()(const double& s, const double& p) const;

 private:
    int m_mom;
    BlattWeisskopf* m_ff;
    double m_precalc;
};

#endif  // SRC_BWWIDTH_H_
