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
    std::unique_ptr<HighFive::DataSet> magnetisation_set;
    std::unique_ptr<HighFive::DataSet> energy_set;
    std::unique_ptr<HighFive::DataSet> rg_magnetisation_set1;
    std::unique_ptr<HighFive::DataSet> rg_magnetisation_set2;
    std::unique_ptr<HighFive::DataSet> rg_magnetisation_set3;
    std::unique_ptr<HighFive::DataSet> rg_energy_set1;
    std::unique_ptr<HighFive::DataSet> rg_energy_set2;
    std::unique_ptr<HighFive::DataSet> rg_energy_set3;
    std::unique_ptr<HighFive::DataSet> cluster_size_set;
    std::string model_type;
    std::filesystem::path project_folder_path;
    int L;
    int dim;
    size_t N;
    double T;
    double beta;
    double J;
    int potts_q = 1;
    // double H;
    // std::vector<double> vec_H;
    // Always records the last X steps
    size_t total_sweeps;
    // bool record_correlation_length = false;
    // bool record_correlation_function = false;
    bool rg_method = false;
};


#endif
