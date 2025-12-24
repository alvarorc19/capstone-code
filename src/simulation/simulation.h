#ifndef SIMULATION_H
#define SIMULATION_H
#include "simulation/simulationParameters.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <hdf5.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

//include highfive and toml here
//for both reading and writing

template <typename T>
class Simulation {
    private:
        SimulationParameters parameters;
        HighFive::File file;
        int time_step;

    public:
        Simulation(): {}
        void parse_parameters();
        void run();
        void do_step();
        void write();
};

#endif
