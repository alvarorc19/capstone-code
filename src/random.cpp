#include "random.h"
#include <random>
#include <iostream>
#include <vector>

const double PI = 3.141592653589793238463;

double rng::random_real_number(std::mt19937 & engine) {
    std::uniform_real_distribution<double> uniform_dist(0,1);
    return uniform_dist(engine);
}
int rng::random_int_number(std::mt19937 & engine, int start, int end){
    std::uniform_int_distribution<int> uniform_dist(start,end);
    return uniform_dist(engine);
}
double rng::random_angle(std::mt19937 & engine){
    std::uniform_real_distribution<double> uniform_dist(0,2 * PI);
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

std::vector<double> rng::random_angle_array(std::mt19937 &engine, int length) {
    std::vector<double> array;
    std::uniform_real_distribution<double> uniform_dist(0,2 * PI);

    for (int i = 0; i< length; i++) {
        array.push_back(uniform_dist(engine));
    }

    return array;
}
