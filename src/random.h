#ifndef RANDOM_H
#define RANDOM_H
#include <iostream>
#include <random>
#include <vector>

namespace rng {
    // Be careful with seed!!!
    inline std::mt19937 engine(0);
    void update_seed(int new_seed);
    int random_int_number(std::mt19937 &engine, int start = 0, int end = 1);
    double random_angle(std::mt19937 &engine);
    double random_real_number(std::mt19937 &engine);
    std::vector<int> random_int_array(std::mt19937 &engine, int length,
                                      int start = 0, int end = 1);
    std::vector<double> random_angle_array(std::mt19937 &engine, int length);
} // namespace rng

#endif
