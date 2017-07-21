#ifndef RNDMI_H
#define RNDMI_H

#include <cstdint>
#include <random>
#include <chrono>  // seed generator

/** @brief Random integers generator (interface) */
class RndmI {
    int32_t m_lo;
    int32_t m_hi;
    static std::default_random_engine rng;
    std::uniform_int_distribution<int32_t> rndm;

 public:
    RndmI(int32_t lo, int32_t hi, int32_t seed = 0) :
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

std::default_random_engine RndmI::rng;

#endif // RNDMI_H
