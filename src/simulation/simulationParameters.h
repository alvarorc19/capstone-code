#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
// Include highfive?
//

//*********************
// lattice_dataset
// model type
// L
// dim
// J
// temp
// potts_q
// time_steps
// recording_steps
// specifying what to record
//   - magnetisation
//   - energy
//   - susceptibility
//   - specific heat
//   - correlation length
//   - correlation function


struct SimulationParameters {
    std::unique_ptr<HighFive::DataSet> lattice_set;
    std::string model_type;
    std::filesystem::path project_folder_path;
    int L;
    int dim;
    double T;
    double beta;
    double J;
    int potts_q;
    int time_steps;
    double H;
    std::vector<double> vec_H;
    // Always records the last X steps
    size_t recording_steps;
    bool record_magnetisation = true;
    bool record_energy = true;
    bool record_susceptibility = false;
    bool record_specific_heat = false;
    bool record_correlation_length = false;
    bool record_correlation_function = false;
};


#endif
