#ifndef SIMULATION_H
#define SIMULATION_H
#include "simulation/simulationParameters.h"
#include "models/modelBase.h"
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <hdf5.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

//include highfive and toml here
//for both reading and writing

// template <typename T>
class Simulation {
    private:
        SimulationParameters parameters;
        HighFive::File file;
        int time_step = 0;
        std::unique_ptr<ModelBase> model;

    public:
        Simulation():parameters{}{}
        void parse_parameters(std::string project_folder);
        void initialise_model();
        void initialise_writing();
        void run();
        void do_metropolis_step();
        void do_cluster_step();
        void write();
};

#endif
