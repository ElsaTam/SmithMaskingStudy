#pragma once

#include <random>
#include <chrono>

struct ScalarRng
{
    std::mt19937_64 rng;
    std::uniform_real_distribution<float> unif;

    inline ScalarRng()
    {
        init(0, 1);
    }

    inline ScalarRng(float min, float max)
    {
        init(min, max);
    }

    inline void init(float min = 0, float max = 1)
    {
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
        rng.seed(ss);
        setRange(min, max);
    }

    inline void setRange(float min, float max)
    {
        unif = std::uniform_real_distribution<float>(min, max);
    }

    // Generate random real in [unif.a(), unif.b()[
    inline float operator() ()
    {
        float f = unif(rng);
        while (f >= unif.b() - m_eps_5)
            f = unif(rng);
        return f;
    }
};


template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}
