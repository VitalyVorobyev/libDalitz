/** Copyright 2017 Vitaly Vorobyev
 ** @file dstrip.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DSTRIP_H_
#define SRC_DSTRIP_H_

///
/// \brief The DStrip class. Defines a strip in the Dalitz phase space
///
class DStrip {
 public:
    DStrip(const double& l, const double& r, const int t) :
        ledge(l), redge(r), type(t) {}
    DStrip(void) : DStrip(0, 0, 0) {}

    double mean(void) const {return 0.5 * (redge + ledge);}
    double sigma(void) const {return 0.5 * (redge - ledge);}

    bool IsWithin(const double& mABsq, const double& mACsq,
                  const double& mBCsq) const {
        switch (type) {
        case 1:  // AB
            if (mABsq > ledge && mABsq < ledge) return true;
            break;
        case 2:  // AC
            if (mACsq > ledge && mACsq < ledge) return true;
            break;
        case 3:  // BC
            if (mBCsq > ledge && mBCsq < ledge) return true;
            break;
        }
        return false;
    }
    double ledge;
    double redge;
    int type;  // 1 -> AB, 2 -> AC, 3 -> BC
};

#endif  // SRC_DSTRIP_H_
