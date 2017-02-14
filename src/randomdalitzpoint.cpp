/** Copyright 2017 Vitaly Vorobyev
 ** @file randomdalitzpoint.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/randomdalitzpoint.h"

#include <iostream>
#include <boost/chrono/chrono.hpp>

typedef std::uniform_real_distribution<double> unif;
typedef std::normal_distribution<double>       gaus;
typedef std::vector<double> vectd;

typedef RandomDalitzPoint RDPoint;

using std::cout;
using std::endl;

std::default_random_engine RandomDalitzPoint::re;
unsigned RandomDalitzPoint::m_seed = 1;

RandomDalitzPoint::RandomDalitzPoint(const double& mmo, const double& mca,
                                     const double& mcb, const double& mcc) :
    DalitzPhaseSpace(mmo, mca, mcb, mcc),
    m_max_tries(1e6),
    mABsqMin(mABsq_min()), mABsqMax(mABsq_max()),
    mACsqMin(mACsq_min()), mACsqMax(mACsq_max()) {
    init();
}

RandomDalitzPoint::RandomDalitzPoint(const DalitzPhaseSpace *phsp) :
    RandomDalitzPoint(phsp->mM(), phsp->mA(), phsp->mB(), phsp->mC()) {}

void RDPoint::SetSeed(const unsigned seed) {
    m_seed = seed;
    re.seed(m_seed);
}

void RDPoint::init(void) {
    const auto now = boost::chrono::system_clock::now().time_since_epoch();
    SetSeed(boost::chrono::duration_cast
            <boost::chrono::milliseconds>(now).count());
}

int RDPoint::GetPoint(double* mABsq, double* mACsq) const {
    unif unifAB(mABsqMin, mABsqMax);
    unif unifAC(mACsqMin, mACsqMax);
    uint64_t ntries = 0;
    do {
        *mABsq = unifAB(re); *mACsq = unifAC(re);
        if (++ntries > m_max_tries) {
            cout << "GetStripPoint: Tries limit exceeded" << endl;
            return -1;
        }
    } while (!IsInPlot(*mABsq, *mACsq));
    return 0;
}

int RDPoint::GetUnconstrainedPoint(double* mABsq, double* mACsq) const {
    unif unifAB(mABsqMin, mABsqMax);
    unif unifAC(mACsqMin, mACsqMax);
    *mABsq = unifAB(re); *mACsq = unifAC(re);
    return 0;
}

int RDPoint::GetPoints(vectd *mABsqV, vectd *mACsqV, const unsigned N) const {
    mABsqV->clear(); mACsqV->clear();
    double mABsq, mACsq;
    for (unsigned i=0; i < N; i++) {
        if (GetPoint(&mABsq, &mACsq)) return -1;
        mABsqV->push_back(mABsq);
        mACsqV->push_back(mACsq);
    }
    return 0;
}

int RDPoint::GetStripPoint(double* mABsq, double* mACsq,
                           const DStrip* shape) const {
    unif unifStrip(shape->mean(), shape->sigma());
    if (GetUnifPoint(mABsq, mACsq, shape, &unifStrip)) return -1;
    return 0;
}

int RDPoint::GetStripPoints(vectd* mABsqV, vectd* mACsqV, const unsigned N,
                            const DStrip* shape) const {
    mABsqV->clear(); mACsqV->clear();
    double mABsq, mACsq;
    unif unifStrip(shape->mean(), shape->sigma());
    for (unsigned i=0; i < N; i++) {
        if (GetUnifPoint(&mABsq, &mACsq, shape, &unifStrip)) return -1;
        mABsqV->push_back(mABsq); mACsqV->push_back(mACsq);
    }
    return 0;
}

int RDPoint::GetStripPoint(double* mABsq, double* mACsq,
                           const DStrip* shape1, const DStrip* shape2) const {
    if (shape1->type == shape2->type) return -5;
    unif unifStrip1(shape1->mean(), shape1->sigma());
    unif unifStrip2(shape2->mean(), shape2->sigma());

    unif& un1 = shape1->type < shape2->type ? unifStrip1 : unifStrip2;
    unif& un2 = shape1->type < shape2->type ? unifStrip2 : unifStrip1;
    const DStrip* sh1 = shape1->type < shape2->type ? shape1 : shape2;
    const DStrip* sh2 = shape1->type < shape2->type ? shape2 : shape1;

    if (GetUnifPoint(mABsq, mACsq, sh1, &un1, sh2, &un2)) return -1;
    return 0;
}

int RDPoint::GetStripPoints(vectd* mABsqV, vectd* mACsqV, const unsigned N,
                            const DStrip* shape1, const DStrip* shape2) const {
    mABsqV->resize(N); mACsqV->resize(N);
    if (shape1->type == shape2->type) return -5;
    double mABsq, mACsq;
    unif unifStrip1(shape1->mean(), shape1->sigma());
    unif unifStrip2(shape2->mean(), shape2->sigma());

    unif& un1 = shape1->type < shape2->type ? unifStrip1 : unifStrip2;
    unif& un2 = shape1->type < shape2->type ? unifStrip2 : unifStrip1;
    const DStrip* sh1 = shape1->type < shape2->type ? shape1 : shape2;
    const DStrip* sh2 = shape1->type < shape2->type ? shape2 : shape1;

    for (unsigned i=0; i < N; i++) {
        if (GetUnifPoint(&mABsq, &mACsq, sh1, &un1, sh2, &un2)) return -1;
        mABsqV->at(i) = mABsq; mACsqV->at(i) = mACsq;
    }
    return 0;
}

int RDPoint::GetGaussPoint(double* mABsq, double* mACsq, double* mBCsq,
                           const DStrip* shape) const {
    gaus GausStrip(shape->mean(), shape->sigma());
    if (GetGausPoint1(mABsq, mACsq, mBCsq, shape, &GausStrip)) return -1;
    return 0;
}

int RDPoint::GetGaussPoints(vectd* mABsqV, vectd* mACsqV, vectd* mBCsqV,
                            const unsigned N, const DStrip* shape) const {
    mABsqV->resize(N); mACsqV->resize(N); mBCsqV->resize(N);
    double mABsq, mACsq, mBCsq;
    gaus GausStrip(shape->mean(), shape->sigma());
    for (unsigned i=0; i < N; i++) {
        if (GetGausPoint1(&mABsq, &mACsq, &mBCsq, shape, &GausStrip)) return -1;
        mABsqV->at(i) = mABsq; mACsqV->at(i) = mACsq; mBCsqV->at(i) = mBCsq;
    }
    return 0;
}

int RDPoint::GetGaussPoint(double* mABsq, double* mACsq, double* mBCsq,
                           const DStrip* shape1, const DStrip* shape2) const {
    if (shape1->type == shape2->type) return -5;
    gaus GausStrip1(shape1->mean(), shape1->sigma());
    gaus GausStrip2(shape2->mean(), shape2->sigma());

    gaus& un1 = shape1->type < shape2->type ? GausStrip1 : GausStrip2;
    gaus& un2 = shape1->type < shape2->type ? GausStrip2 : GausStrip1;
    const DStrip* sh1 = shape1->type < shape2->type ? shape1 : shape2;
    const DStrip* sh2 = shape1->type < shape2->type ? shape2 : shape1;

    if (GetGausPoint2(mABsq, mACsq, mBCsq, sh1, &un1, sh2, &un2)) return -1;
    return 0;
}

int RDPoint::GetGaussPoints(vectd* mABsqV, vectd* mACsqV, vectd* mBCsqV,
                            const unsigned N,
                            const DStrip* shape1, const DStrip* shape2) const {
    mABsqV->resize(N); mACsqV->resize(N); mBCsqV->resize(N);
    if (shape1->type == shape2->type) return -5;
    double mABsq, mACsq, mBCsq;
    gaus GausStrip1(shape1->mean(), shape1->sigma());
    gaus GausStrip2(shape2->mean(), shape2->sigma());

    gaus& un1 = shape1->type < shape2->type ? GausStrip1 : GausStrip2;
    gaus& un2 = shape1->type < shape2->type ? GausStrip2 : GausStrip1;
    const DStrip* sh1 = shape1->type < shape2->type ? shape1 : shape2;
    const DStrip* sh2 = shape1->type < shape2->type ? shape2 : shape1;

    for (unsigned i=0; i < N; i++) {
        if (GetGausPoint2(&mABsq, &mACsq, &mBCsq, sh1, &un1, sh2, &un2))
            return -1;
        mABsqV->at(i) = mABsq; mACsqV->at(i) = mACsq; mBCsqV->at(i) = mBCsq;
    }
    return 0;
}

// Private methods
int RDPoint::GetUnifPoint(double* mABsq, double* mACsq,
                          const DStrip* shape, unif* dist) const {
    uint64_t ntries = 0;
    unif unifAB(mABsqMin, mABsqMax);
    unif unifAC(mACsqMin, mACsqMax);
    double mBCsq;
    switch (shape->type) {
    case 1:  // AB
        do {
            if (++ntries > m_max_tries) return -1;
            *mABsq = (*dist)(re);
            *mACsq = unifAC(re);
            mBCsq = GetmBCsq(*mABsq, *mACsq);
        } while (!shape->IsWithin(*mABsq, *mACsq, mBCsq) ||
                 !IsInPlot(*mABsq, *mACsq));
        break;
    case 2:  // AC
        do {
            if (++ntries > m_max_tries) return -2;
            *mABsq = unifAB(re);
            *mACsq = (*dist)(re);
            mBCsq = GetmBCsq(*mABsq, *mACsq);
        } while (!shape->IsWithin(*mABsq, *mACsq, mBCsq) ||
                 !IsInPlot(*mABsq, *mACsq));
        break;
    case 3:  // AC
        do {
            if (++ntries > m_max_tries) return -3;
            *mABsq = unifAB(re);
            mBCsq = (*dist)(re);
            *mACsq = GetmBCsq(*mABsq, mBCsq);
        } while (!shape->IsWithin(*mABsq, *mACsq, mBCsq) ||
                 !IsInPlot(*mABsq, *mACsq));
        break;
    default:
        return -4;
    }
    return 0;
}

int RDPoint::GetGausPoint1(double* mABsq, double* mACsq, double* mBCsq,
                           const DStrip* shape, gaus* dist) const {
    uint64_t ntries = 0;
    unif unifAB(mABsqMin, mABsqMax);
    unif unifAC(mACsqMin, mACsqMax);
    switch (shape->type) {
    case 1:  // AB
        do {
            if (++ntries > m_max_tries) return -1;
            *mABsq = (*dist)(re);
            *mACsq = unifAC(re);
            *mBCsq = GetmBCsq(*mABsq, *mACsq);
        } while (!IsInPlot(*mABsq, *mACsq));
        break;
    case 2:  // AC
        do {
            if (++ntries > m_max_tries) return -2;
            *mABsq = unifAB(re);
            *mACsq = (*dist)(re);
            *mBCsq = GetmBCsq(*mABsq, *mACsq);
        } while (!IsInPlot(*mABsq, *mACsq));
        break;
    case 3:  // AC
        do {
            if (++ntries > m_max_tries) return -3;
            *mABsq = unifAB(re);
            *mBCsq = (*dist)(re);
            *mACsq = GetmBCsq(*mABsq, *mBCsq);
        } while (!IsInPlot(*mABsq, *mACsq));
        break;
    default:
        return -4;
    }
    return 0;
}

int RDPoint::GetUnifPoint(double* mABsq, double* mACsq,
                          const DStrip* shape1, unif* dist1,
                          const DStrip* shape2, unif* dist2) const {
    if (shape1->type >= shape2->type) return -5;
    uint64_t ntries = 0;
    double mBCsq;
    switch (shape1->type) {
    case 1:  // AB
        if (shape2->type == 2) {  // AC
            do {
                if (++ntries > m_max_tries) return -1;
                *mABsq = (*dist1)(re);
                *mACsq = (*dist2)(re);
                mBCsq = GetmBCsq(*mABsq, *mACsq);
            } while (!shape1->IsWithin(*mABsq, *mACsq, mBCsq) ||
                     !shape2->IsWithin(*mABsq, *mACsq, mBCsq) ||
                     !IsInPlot(*mABsq, *mACsq));
        } else {  // BC
            do {
                if (++ntries > m_max_tries) return -2;
                *mABsq = (*dist1)(re);
                mBCsq = (*dist2)(re);
                *mACsq = GetmBCsq(*mABsq, mBCsq);
            } while (!shape1->IsWithin(*mABsq, *mACsq, mBCsq) ||
                     !shape2->IsWithin(*mABsq, *mACsq, mBCsq) ||
                     !IsInPlot(*mABsq, *mACsq));
        }
        break;
    case 2:  // AC -> shape2.type == 3
        do {
            if (++ntries > m_max_tries) return -3;
            *mACsq = (*dist1)(re);
            mBCsq = (*dist2)(re);
            *mABsq = GetmBCsq(*mACsq, mBCsq);
        } while (!shape1->IsWithin(*mABsq, *mACsq, mBCsq) ||
                 !shape2->IsWithin(*mABsq, *mACsq, mBCsq) ||
                 !IsInPlot(*mABsq, *mACsq));
    }
    return 0;
}

int RDPoint::GetGausPoint2(double* mABsq, double* mACsq, double* mBCsq,
                           const DStrip* shape1, gaus* dist1,
                           const DStrip* shape2, gaus* dist2) const {
    if (shape1->type >= shape2->type) return -1;
    uint64_t ntries = 0;
    switch (shape1->type) {
    case 1:  // AB
        if (shape2->type == 2) {  // AC
            do {
                if (++ntries > m_max_tries) return -1;
                *mABsq = (*dist1)(re);
                *mACsq = (*dist2)(re);
                *mBCsq = GetmBCsq(*mABsq, *mACsq);
            } while (!IsInPlot(*mABsq, *mACsq));
        } else {  // BC
            do {
                if (++ntries > m_max_tries) return -2;
                *mABsq = (*dist1)(re);
                *mBCsq = (*dist2)(re);
                *mACsq = GetmBCsq(*mABsq, *mBCsq);
            } while (!IsInPlot(*mABsq, *mACsq));
        }
        break;
    case 2:  // AC -> shape2.type == 3
        do {
            if (++ntries > m_max_tries) return -3;
            *mACsq = (*dist1)(re);
            *mBCsq = (*dist2)(re);
            *mABsq = GetmBCsq(*mACsq, *mBCsq);
        } while (!IsInPlot(*mABsq, *mACsq));
    }
    return 0;
}
