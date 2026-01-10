#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H
#include <iostream>
#include <string>
#include <vector>
#include <hdf5.h>
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
    HighFive::DataSet lattice_set;
    std::string model_type;
    int L;
    int dim;
    double T;
    double J;
    int potts_q;
    int time_steps;
    // Always records the last X steps
    int recording_steps;
    bool record_magnetisation = true;
    bool record_energy = true;
    bool record_susceptibility = false;
    bool record_specific_heat = false;
    bool record_correlation_length = false;
    bool record_correlation_function = false;
};


#endif
