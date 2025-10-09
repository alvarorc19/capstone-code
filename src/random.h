#ifndef RANDOM_
#define RANDOM_
#include <iostream>
#include <random>
#include <vector>

namespace rng {
const int seed = 42;
inline std::mt19937 engine(seed);
double random_real_number(std::mt19937 &engine);
std::vector<int> random_int_array(std::mt19937 &engine, int length,
                                  int start = 0, int end = 1);
} // namespace rng

#endif
