#include <random>
#include <iostream>
#include <vector>
#include "random.h"

double rng::random_real_number(std::mt19937 & engine) {
    std::uniform_real_distribution<double> uniform_dist(0,1);
    return uniform_dist(engine);
}

std::vector<int> rng::random_int_array(std::mt19937 & engine, int length, int start, int end) {
    std::vector<int> int_array;
    std::uniform_int_distribution<int> uniform_dist(start,end);

    for (int i = 0; i<length; ++i) {
        int_array.push_back(uniform_dist(engine));
    }

    return int_array;
}
