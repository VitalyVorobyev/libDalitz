#ifndef DSTRIP_H
#define DSTRIP_H

typedef const double cdouble;

///
/// \brief The DStrip class. Defines a strip in the Dalitz phase space
///
class DStrip{
public:
  DStrip(cdouble& l, cdouble& r, const int t): ledge(l),redge(r),type(t) {}
  DStrip(void): DStrip(0,0,0) {}

  double mean(void)  const {return 0.5*(redge+ledge);}
  double sigma(void) const {return 0.5*(redge-ledge);}

  bool IsWithin(cdouble& mABsq, cdouble& mACsq, cdouble& mBCsq) const{
    switch(type){
    case 1:// AB
      if(mABsq>ledge && mABsq<ledge) return true;
      break;
    case 2:// AC
      if(mACsq>ledge && mACsq<ledge) return true;
      break;
    case 3:// BC
      if(mBCsq>ledge && mBCsq<ledge) return true;
      break;
    }
    return false;
  }
  double ledge;
  double redge;
  int type;// 1 -> AB, 2 -> AC, 3 -> BC
};

#endif // DSTRIP_H
