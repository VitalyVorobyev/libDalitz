#ifndef RNDMD_H
#define RNDMD_H

#include <cstdint>
#include <random>
#include <chrono>  // seed generator

/** @brief Random double generator (interface) */
class RndmD {
    double m_lo;
    double m_hi;
    static std::default_random_engine rng;
    std::uniform_real_distribution<double> rndm;

 public:
    RndmD(double lo, double hi, int32_t seed = 0) :
        m_lo(lo), m_hi(hi), rndm(m_lo, m_hi) {
        if (!seed)
            rng.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        else rng.seed(seed);
    }

    auto operator()() {return rndm(rng);}
    auto nonzero() {
        auto var = 0;
        while (!(var = rndm(rng))) {}
        return var;}
};

std::default_random_engine RndmD::rng;

#endif // RNDMD_H
