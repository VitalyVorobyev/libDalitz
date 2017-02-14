/** Copyright 2017 Vitaly Vorobyev
 ** @file randomdalitzpoint.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_RANDOMDALITZPOINT_H_
#define SRC_RANDOMDALITZPOINT_H_

#include <random>  // std::default_random_engine
#include <vector>
#include <cstdint>

#include "./dalitzphasespace.h"
#include "./dstrip.h"

///
/// \brief The RandomDalitzPoint class
///
class RandomDalitzPoint : public DalitzPhaseSpace {
 public:
    explicit RandomDalitzPoint(const DalitzPhaseSpace* phsp);
    RandomDalitzPoint(const double& mmo, const double& mca,
                      const double& mcb, const double& mcc);
    int GetPoint(double* mABsq, double* mACsq) const;
    int GetUnconstrainedPoint(double* mABsq, double* mACsq) const;
    int GetStripPoint(double* mABsq, double* mACsq, const DStrip* shape) const;
    int GetGaussPoint(double* mABsq, double* mACsq, double* mBCsq,
                      const DStrip* shape) const;
    int GetStripPoint(double* mABsq, double* mACsq, const DStrip* shape1,
                      const DStrip* shape2) const;
    int GetGaussPoint(double* mABsq, double* mACsq, double* mBCsq,
                      const DStrip* shape1, const DStrip* shape2) const;

    int GetPoints(std::vector<double>* mABsqV, std::vector<double>* mACsqV,
                  const unsigned N) const;
    int GetStripPoints(std::vector<double>* mABsqV,
                       std::vector<double>* mACsqV,
                       const unsigned N, const DStrip* shape) const;
    int GetGaussPoints(std::vector<double>* mABsqV,
                       std::vector<double>* mACsqV,
                       std::vector<double>* mBCsqV,
                       const unsigned N, const DStrip* shape) const;
    int GetStripPoints(std::vector<double>* mABsqV,
                       std::vector<double>* mACsqV,
                       const unsigned N, const DStrip* shape1,
                       const DStrip* shape2) const;
    int GetGaussPoints(std::vector<double>* mABsqV,
                       std::vector<double>* mACsqV,
                       std::vector<double>* mBCsqV,
                       const unsigned N, const DStrip* shape1,
                       const DStrip* shape2) const;
    // * Static methods * //
    static void SetSeed(const unsigned seed);
    static unsigned GetSeed(void) {return m_seed;}

 private:
    void init();

    int GetUnifPoint(double* mABsq, double* mACsq, const DStrip* shape,
                     std::uniform_real_distribution<double>* dist) const;
    int GetGausPoint1(double* mABsq, double* mACsq, double* mBCsq,
                      const DStrip* shape,
                      std::normal_distribution<double>* dist) const;
    int GetUnifPoint(double* mABsq, double* mACsq, const DStrip* shape1,
                     std::uniform_real_distribution<double>* dist1,
                     const DStrip* shape2,
                     std::uniform_real_distribution<double>* dist2) const;
    int GetGausPoint2(double* mABsq, double* mACsq, double* mBCsq,
                      const DStrip* shape1,
                      std::normal_distribution<double>* dist1,
                      const DStrip* shape2,
                      std::normal_distribution<double>* dist2) const;
    uint64_t m_max_tries;
    const double mABsqMin;
    const double mABsqMax;
    const double mACsqMin;
    const double mACsqMax;

    static std::default_random_engine re;
    static unsigned m_seed;
};

#endif  // SRC_RANDOMDALITZPOINT_H_
