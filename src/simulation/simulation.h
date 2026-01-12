#ifndef SIMULATION_H
#define SIMULATION_H
#include "simulation/simulationParameters.h"
#include "simulation/observables.h"
#include "models/modelBase.h"
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <toml++/toml.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

//include highfive and toml here
//for both reading and writing

// template <typename T>
class Simulation {
    private:
        SimulationParameters parameters;
        Observables observables;
        std::unique_ptr<HighFive::File> file;
        int time_step = 0;
        std::unique_ptr<ModelBase> model;

        ivec extract_int_vector_toml(const toml::array& array);
        dvec extract_double_vector_toml(const toml::array& array);

    public:
        Simulation() =default;
        void parse_parameters(std::string project_folder, std::string model_type);
        void initialise_model();
        void initialise_writing();
        void run();
        void do_metropolis_step();
        void do_cluster_step();
        void write_lattice(int time);
        void write_observables();
        void update_observables();
};

#endif
