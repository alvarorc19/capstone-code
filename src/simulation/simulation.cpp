#include "simulation/simulation.h"
#include "simulation/simulationParameters.h"
#include "simulation/simulationParameters.h"
#include "models/ising.h"
#include "models/potts.h"
#include "models/xy.h"
#include "random.h"

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <omp.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <toml++/toml.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

void Simulation::parse_parameters(std::filesystem::path project_folder_path, std::string model_type) {
    toml::table tbl;
    // std::filesystem::path project_folder_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
    // careful with this, has not been tested
    try
    {
        std::filesystem::path config_path = project_folder_path / std::filesystem::path("config.toml");
        tbl = toml::parse_file(config_path.string());
        // std::cout << tbl << "\n";
    }
    catch (const toml::parse_error& err)
    {
        std::cout << "Make sure the parameters are in a file named config.toml" << "\n";
        std::cerr << "Parsing failed:\n" << err << "\n";
        // return 2;
    }

    // Assign the seed
    int seed = static_cast<int>(tbl["seed"].value_or<int64_t>(0));
    rng::update_seed(seed);

    // Assign all the parameters to the struct
    parameters.project_folder_path = project_folder_path;
    parameters.dim = static_cast<int>(tbl["physical_settings"]["dimension"].value_or<int64_t>(0));
    parameters.L = static_cast<int>(tbl["physical_settings"]["L"].value_or<int64_t>(0));
    parameters.N = std::pow(parameters.L, parameters.dim);
    parameters.T = tbl["physical_settings"]["temperature"].value_or<double>(0.0);
    parameters.beta = 1 / parameters.T;
    parameters.J = tbl["physical_settings"]["J"].value_or<double>(0.0);
    // parameters.potts_q = static_cast<int>(tbl["physical_settings"]["potts_q"].value_or<int64_t>(0));
    // parameters.H = tbl["physical_settings"]["H"].value_or<double>(0.0);
    //
    // parameters.vec_H = extract_double_vector_toml(*tbl["physical_settings"]["vec_H"].as_array());

    parameters.model_type = model_type;
    parameters.total_sweeps = static_cast<size_t>(tbl["simulation_settings"]["total_sweeps"].value_or<int64_t>(0));
    // parameters.record_correlation_length = tbl["simulation_settings"]["record_correlation_length"].value_or<bool>(false);
    // parameters.record_correlation_function = tbl["simulation_settings"]["record_correlation_function"].value_or<bool>(false);
    parameters.rg_method = tbl["simulation_settings"]["rg_method"].value_or<bool>(false);
}

void Simulation::run(){
    if (parameters.rg_method and parameters.dim == 2) {
        rg_run_2d();
    }
    else if (parameters.rg_method and parameters.dim == 3) {
        rg_run_3d();
    }

    else {
        normal_run();
    }

}

void Simulation::normal_run() {

    int num_loops = parameters.total_sweeps >> 9;
    int modulo = parameters.total_sweeps % 512;

    if (modulo == 0){
        num_loops -= 1;
        modulo += 512;
    }

    visited.assign(parameters.N, 0);
    cluster_stack.reserve(parameters.N);
    cluster_size.reserve(parameters.N);
    
    // Start timer
    const auto start{std::chrono::steady_clock::now()};

    do_cluster_sweep();
    update_observables(0, cluster_size);
    time_step++;

    // End timer
    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};

    #pragma omp critical(io)
    {
        std::cout << "One sweep took " << elapsed_seconds.count() << " s\n";
        std::cout << "Expected execution time (not accounting for writing) =  " << 
            elapsed_seconds.count() * parameters.total_sweeps << " s " << std::endl;
    }

    for(int j = 1; j < 512; j++){
        do_cluster_sweep();
        update_observables(j, cluster_size);
        time_step++;
    }

    #pragma omp critical (HDF5)
    {
        write_observables_after_loop(0);
    }

    // End timer
    const auto finish512{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds512{finish512 - start};

    #pragma omp critical(io)
    {
        std::cout << "One loop of 512 sweeps took " << elapsed_seconds512.count() <<"s" << "\n";
        std::cout << "Need to do " << num_loops << " loops" << std::endl;
    }

    for (int i = 1; i < num_loops; i++){
        for(int j = 0; j < 512; j++){
            do_cluster_sweep();
            update_observables(j, cluster_size);
            time_step++;
        }

        #pragma omp critical (HDF5)
        {
            write_observables_after_loop(i << 9);
        }
    }

    dvec final_x_magnetisation;
    dvec final_y_magnetisation;
    dvec final_energy;
    dvec final_average_cluster_size;
    final_x_magnetisation.reserve(modulo);
    final_y_magnetisation.reserve(modulo);
    final_energy.reserve(modulo);
    final_average_cluster_size.reserve(modulo);

    for (int i = 0; i < modulo - 1; i++){
        do_cluster_sweep();
        final_energy.emplace_back(model->compute_total_energy());
        final_x_magnetisation.emplace_back(model->compute_spin_magnetic_term(0));
        final_y_magnetisation.emplace_back(model->compute_spin_magnetic_term(1));
        final_average_cluster_size.emplace_back(compute_average_cluster_size(cluster_size));
        time_step++;
    }


    // Record last 100 steps of the lattice
    cluster_size.clear();
    spins_flipped = 0;
    int flip;
    for (int i = 0; i < parameters.N - 100 ; i++){
        flip = spins_flipped; 
        do_cluster_step();
        cluster_size.emplace_back(spins_flipped - flip);
    }

    for (int i = 0; i < 100; i++){
        flip = spins_flipped; 
        do_cluster_step();
        cluster_size.emplace_back(spins_flipped - flip);
        #pragma omp critical(HDF5)
        {
            write_lattice(i);
        }
    }
    
    final_energy.emplace_back(model->compute_total_energy());
    final_x_magnetisation.emplace_back(model->compute_spin_magnetic_term(0));
    final_y_magnetisation.emplace_back(model->compute_spin_magnetic_term(1));
    final_average_cluster_size.emplace_back(compute_average_cluster_size(cluster_size));
    time_step++;

    #pragma omp critical (HDF5)
    {
        write_observables_final(num_loops << 9, final_energy, final_x_magnetisation,final_y_magnetisation, final_average_cluster_size);
        file->flush();
        close_writing();
    }

    #pragma omp critical(io)
    {
        std::cout << "Finished run for " << parameters.project_folder_path.filename() << "\n";
    }
}

void Simulation::rg_run_2d() {

    int num_loops = parameters.total_sweeps >> 9;
    int modulo = parameters.total_sweeps % 512;

    int b1 = 2;
    int b2 = 4;
    int b3 = 16;
    ivec neighbours_table1 = model->calculate_reduced_neighbours_table(parameters.L, parameters.dim, b1);
    ivec neighbours_table2 = model->calculate_reduced_neighbours_table(parameters.L, parameters.dim, b2);
    ivec neighbours_table3 = model->calculate_reduced_neighbours_table(parameters.L, parameters.dim, b3);

    if (modulo == 0){
        num_loops -= 1;
        modulo += 512;
    }

    visited.assign(parameters.N, 0);
    cluster_stack.reserve(parameters.N);
    cluster_size.reserve(parameters.N);
    
    // Start timer
    const auto start{std::chrono::steady_clock::now()};

    do_cluster_sweep();
    update_observables(0, cluster_size);
    update_rg_observables(0, b1, b2,b3, neighbours_table1, neighbours_table2, neighbours_table3);
    time_step++;

    // End timer
    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};

    #pragma omp critical(io)
    {
        std::cout << "One sweep took " << elapsed_seconds.count() << " s\n";
        std::cout << "Expected execution time (not accounting for writing) =  " << 
            elapsed_seconds.count() * parameters.total_sweeps << " s " << std::endl;
    }

    for(int j = 1; j < 512; j++){
        do_cluster_sweep();
        update_observables(j, cluster_size);
        update_rg_observables(j, b1, b2,b3, neighbours_table1, neighbours_table2, neighbours_table3);
        time_step++;
    }

    #pragma omp critical (HDF5)
    {
        write_observables_after_loop(0);
        write_rg_observables_after_loop(0);
    }

    // End timer
    const auto finish512{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds512{finish512 - start};

    #pragma omp critical(io)
    {
        std::cout << "One loop of 512 sweeps took " << elapsed_seconds512.count() <<"s" << "\n";
        std::cout << "Need to do " << num_loops << " loops" << std::endl;
    }

    for (int i = 1; i < num_loops; i++){
        for(int j = 0; j < 512; j++){
            do_cluster_sweep();
            update_observables(j, cluster_size);
            update_rg_observables(j, b1, b2,b3, neighbours_table1, neighbours_table2, neighbours_table3);
            time_step++;
        }

        #pragma omp critical (HDF5)
        {
            write_observables_after_loop(i << 9);
            write_rg_observables_after_loop(i << 9);
        }
    }

    dvec final_x_magnetisation;
    dvec final_y_magnetisation;
    dvec final_energy;
    dvec final_average_cluster_size;
    dvec final_rg_x_magnetisation1;
    dvec final_rg_y_magnetisation1;
    dvec final_rg_x_magnetisation2;
    dvec final_rg_y_magnetisation2;
    dvec final_rg_x_magnetisation3;
    dvec final_rg_y_magnetisation3;
    dvec final_rg_energy1;
    dvec final_rg_energy2;
    dvec final_rg_energy3;
    final_x_magnetisation.reserve(modulo);
    final_y_magnetisation.reserve(modulo);
    final_energy.reserve(modulo);
    final_average_cluster_size.reserve(modulo);
    final_rg_x_magnetisation1.reserve(modulo);
    final_rg_y_magnetisation1.reserve(modulo);
    final_rg_x_magnetisation2.reserve(modulo);
    final_rg_y_magnetisation2.reserve(modulo);
    final_rg_x_magnetisation3.reserve(modulo);
    final_rg_y_magnetisation3.reserve(modulo);
    final_rg_energy1.reserve(modulo);
    final_rg_energy2.reserve(modulo);
    final_rg_energy3.reserve(modulo);

    for (int i = 0; i < modulo - 1; i++){
        do_cluster_sweep();
        final_energy.emplace_back(model->compute_total_energy());
        final_x_magnetisation.emplace_back(model->compute_spin_magnetic_term(0));
        final_y_magnetisation.emplace_back(model->compute_spin_magnetic_term(1));
        final_average_cluster_size.emplace_back(compute_average_cluster_size(cluster_size));

        // RG 1
        model->compute_reduced_lattice(b1, parameters.N, parameters.L);
        final_rg_x_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(0));
        final_rg_y_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(1));
        final_rg_energy1.emplace_back(model->compute_rg_energy(b1, parameters.N, neighbours_table1, parameters.dim));
        // RG 2
        model->compute_reduced_lattice(b2, parameters.N, parameters.L);
        final_rg_x_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(0));
        final_rg_y_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(1));
        final_rg_energy2.emplace_back(model->compute_rg_energy(b2, parameters.N, neighbours_table2, parameters.dim));
        // RG 3
        model->compute_reduced_lattice(b3, parameters.N, parameters.L);
        final_rg_x_magnetisation3.emplace_back(model->compute_rg_spin_magnetic_term(0));
        final_rg_y_magnetisation3.emplace_back(model->compute_rg_spin_magnetic_term(1));
        final_rg_energy3.emplace_back(model->compute_rg_energy(b3, parameters.N, neighbours_table3, parameters.dim));
        time_step++;
    }


    // Record last 100 steps of the lattice
    cluster_size.clear();
    spins_flipped = 0;
    int flip;
    for (int i = 0; i < parameters.N - 100 ; i++){
        flip = spins_flipped; 
        do_cluster_step();
        cluster_size.emplace_back(spins_flipped - flip);
    }

    for (int i = 0; i < 100; i++){
        flip = spins_flipped; 
        do_cluster_step();
        cluster_size.emplace_back(spins_flipped - flip);
        #pragma omp critical(HDF5)
        {
            write_lattice(i);
        }
    }
    
    final_energy.emplace_back(model->compute_total_energy());
    final_x_magnetisation.emplace_back(model->compute_spin_magnetic_term(0));
    final_y_magnetisation.emplace_back(model->compute_spin_magnetic_term(1));
    final_average_cluster_size.emplace_back(compute_average_cluster_size(cluster_size));
    // RG 1
    model->compute_reduced_lattice(b1, parameters.N, parameters.L);
    final_rg_x_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(0));
    final_rg_y_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(1));
    final_rg_energy1.emplace_back(model->compute_rg_energy(b1, parameters.N, neighbours_table1, parameters.dim));
    // RG 2
    model->compute_reduced_lattice(b2, parameters.N, parameters.L);
    final_rg_x_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(0));
    final_rg_y_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(1));
    final_rg_energy2.emplace_back(model->compute_rg_energy(b2, parameters.N, neighbours_table2, parameters.dim));
    // RG 3
    model->compute_reduced_lattice(b3, parameters.N, parameters.L);
    final_rg_x_magnetisation3.emplace_back(model->compute_rg_spin_magnetic_term(0));
    final_rg_y_magnetisation3.emplace_back(model->compute_rg_spin_magnetic_term(1));
    final_rg_energy3.emplace_back(model->compute_rg_energy(b3, parameters.N, neighbours_table3, parameters.dim));
    time_step++;

    #pragma omp critical (HDF5)
    {
        write_observables_final(num_loops << 9, final_energy, final_x_magnetisation,final_y_magnetisation, final_average_cluster_size);
        write_rg_observables_final(num_loops << 9, final_rg_energy1, final_rg_energy2,final_rg_energy3, final_rg_x_magnetisation1, final_rg_x_magnetisation2,final_rg_x_magnetisation3, final_rg_y_magnetisation1, final_rg_y_magnetisation2, final_rg_y_magnetisation3);
        file->flush();
        close_writing();
    }

    #pragma omp critical(io)
    {
        std::cout << "Finished run for " << parameters.project_folder_path.filename() << "\n";
    }
}


void Simulation::initialise_model() {
    if (parameters.model_type == "ising"){
        model = std::make_unique<IsingModel>(
            parameters.T,
            parameters.J,
            // parameters.H,
            parameters.dim,
            parameters.L
        );
    }
    else if (parameters.model_type == "potts"){
        model = std::make_unique<PottsModel>(
            parameters.T,
            parameters.J,
            // parameters.vec_H,
            // parameters.H,
            parameters.dim,
            parameters.L,
            parameters.potts_q
        );
    }
    else if (parameters.model_type == "xy"){
        model = std::make_unique<XYModel>(
            parameters.T,
            parameters.J,
            // parameters.vec_H,
            // parameters.H,
            parameters.dim,
            parameters.L
        );
    }
    else {
        throw std::invalid_argument("Invalid model_type: '" + parameters.model_type + "'. Must be 'ising', 'potts', or 'xy'.");
    }
}

void Simulation::initialise_writing() {
    std::filesystem::path filename = parameters.project_folder_path / std::filesystem::path("results.h5");
    this->file = std::make_unique<HighFive::File>(filename.string(), HighFive::File::Truncate);
    // Define chunks, how data is stored
    std::vector<hsize_t> chunk_dims = {30, parameters.N};
    
    if (parameters.dim == 2){
        // Lattice dataset, stores only last sweep I would do, last N iterations.
        HighFive::DataSpace lattice_space({100, parameters.N});

        HighFive::DataSetCreateProps lattice_props;
        lattice_props.add(HighFive::Chunking(chunk_dims));

        if (parameters.model_type == "xy") {
            parameters.lattice_set = std::make_unique<HighFive::DataSet>(
                file->createDataSet<double>(
                    "lattice",
                    lattice_space,
                    lattice_props
                )
            );
        }
        else if (parameters.model_type == "ising" or parameters.model_type == "potts") {
            parameters.lattice_set = std::make_unique<HighFive::DataSet>(
                file->createDataSet<int>(
                    "lattice",
                    lattice_space,
                    lattice_props
                )
            );
        }
    }

    // Create observables datasets
    HighFive::DataSpace magnetisation_space({parameters.total_sweeps, 2});
    chunk_dims = {512, 2};
    HighFive::DataSetCreateProps magnetisation_props;
    magnetisation_props.add(HighFive::Chunking(chunk_dims));
    parameters.magnetisation_set = std::make_unique<HighFive::DataSet>(
        file->createDataSet<double>(
            "observables/magnetisation",
            magnetisation_space,
            magnetisation_props
        )
    );

    observables.x_magnetisation.resize(512);
    observables.y_magnetisation.resize(512);

    HighFive::DataSpace energy_space({parameters.total_sweeps});
    chunk_dims = {512};
    HighFive::DataSetCreateProps energy_props;
    energy_props.add(HighFive::Chunking(chunk_dims));
    parameters.energy_set = std::make_unique<HighFive::DataSet>(
        file->createDataSet<double>(
            "observables/energy",
            energy_space,
            energy_props
        )
    );

    observables.energy_array.resize(512);

    HighFive::DataSpace cluster_size_space({parameters.total_sweeps});
    chunk_dims = {512};
    HighFive::DataSetCreateProps cluster_size_props;
    cluster_size_props.add(HighFive::Chunking(chunk_dims));
    parameters.cluster_size_set = std::make_unique<HighFive::DataSet>(
        file->createDataSet<double>(
            "observables/average_cluster_size",
            cluster_size_space,
            cluster_size_props
        )
    );

    observables.average_cluster_size.resize(512);

    // TODO: Modularise the stuff to only have the blist (maybe have in
    //      TOML the list?
    if (parameters.rg_method) {
        HighFive::DataSpace rg_magnetisation_space1({parameters.total_sweeps, 2});
        chunk_dims = {512, 2};
        HighFive::DataSetCreateProps rg_magnetisation_props1;
        rg_magnetisation_props1.add(HighFive::Chunking(chunk_dims));
        parameters.rg_magnetisation_set1 = std::make_unique<HighFive::DataSet>(
            file->createDataSet<double>(
                "renormalisation/magnetisation_b_2",
                rg_magnetisation_space1,
                rg_magnetisation_props1
            )
        );

        observables.rg_x_magnetisation1.resize(512);
        observables.rg_y_magnetisation1.resize(512);

        HighFive::DataSpace rg_magnetisation_space2({parameters.total_sweeps, 2});
        chunk_dims = {512, 2};
        HighFive::DataSetCreateProps rg_magnetisation_props2;
        rg_magnetisation_props2.add(HighFive::Chunking(chunk_dims));
        parameters.rg_magnetisation_set2 = std::make_unique<HighFive::DataSet>(
            file->createDataSet<double>(
                "renormalisation/magnetisation_b_4",
                rg_magnetisation_space2,
                rg_magnetisation_props2
            )
        );

        observables.rg_x_magnetisation2.resize(512);
        observables.rg_y_magnetisation2.resize(512);


        HighFive::DataSpace rg_energy_space1({parameters.total_sweeps});
        chunk_dims = {512};
        HighFive::DataSetCreateProps rg_energy_props1;
        rg_energy_props1.add(HighFive::Chunking(chunk_dims));
        parameters.rg_energy_set1 = std::make_unique<HighFive::DataSet>(
            file->createDataSet<double>(
                "renormalisation/energy_b_2",
                rg_energy_space1,
                rg_energy_props1
            )
        );

        observables.rg_energy1.resize(512);

        HighFive::DataSpace rg_energy_space2({parameters.total_sweeps});
        chunk_dims = {512};
        HighFive::DataSetCreateProps rg_energy_props2;
        rg_energy_props2.add(HighFive::Chunking(chunk_dims));
        parameters.rg_energy_set2 = std::make_unique<HighFive::DataSet>(
            file->createDataSet<double>(
                "renormalisation/energy_b_4",
                rg_energy_space2,
                rg_energy_props2
            )
        );

        observables.rg_energy2.resize(512);

        if (parameters.dim == 2){
            HighFive::DataSpace rg_magnetisation_space3({parameters.total_sweeps, 2});
            chunk_dims = {512, 2};
            HighFive::DataSetCreateProps rg_magnetisation_props3;
            rg_magnetisation_props3.add(HighFive::Chunking(chunk_dims));
            parameters.rg_magnetisation_set3 = std::make_unique<HighFive::DataSet>(
                file->createDataSet<double>(
                    "renormalisation/magnetisation_b_16",
                    rg_magnetisation_space3,
                    rg_magnetisation_props3
                )
            );

            observables.rg_x_magnetisation3.resize(512);
            observables.rg_y_magnetisation3.resize(512);

            HighFive::DataSpace rg_energy_space3({parameters.total_sweeps});
            chunk_dims = {512};
            HighFive::DataSetCreateProps rg_energy_props3;
            rg_energy_props3.add(HighFive::Chunking(chunk_dims));
            parameters.rg_energy_set3 = std::make_unique<HighFive::DataSet>(
                file->createDataSet<double>(
                    "renormalisation/energy_b_16",
                    rg_energy_space3,
                    rg_energy_props3
                )
            );

            observables.rg_energy3.resize(512);
        }
    }

}
void Simulation::close_writing() {
    if (parameters.dim == 2){
        parameters.lattice_set.reset();
    }

    parameters.magnetisation_set.reset();
    parameters.energy_set.reset();
    parameters.cluster_size_set.reset();

    if (parameters.rg_method){
        parameters.rg_magnetisation_set1.reset();
        parameters.rg_magnetisation_set2.reset();
        parameters.rg_magnetisation_set3.reset();

        parameters.rg_energy_set1.reset();
        parameters.rg_energy_set2.reset();
        parameters.rg_energy_set3.reset();
    }

    if (file) {
        file->flush();
        file.reset();
    }
}

void Simulation::rg_run_3d() {

    int num_loops = parameters.total_sweeps >> 9;
    int modulo = parameters.total_sweeps % 512;

    int b1 = 2;
    int b2 = 4;
    ivec neighbours_table1 = model->calculate_reduced_neighbours_table(parameters.L, parameters.dim, b1);
    ivec neighbours_table2 = model->calculate_reduced_neighbours_table(parameters.L, parameters.dim, b2);

    if (modulo == 0){
        num_loops -= 1;
        modulo += 512;
    }

    visited.assign(parameters.N, 0);
    cluster_stack.reserve(parameters.N);
    cluster_size.reserve(parameters.N);
    
    // Start timer
    const auto start{std::chrono::steady_clock::now()};

    do_cluster_sweep();
    update_observables(0, cluster_size);
    update_rg_observables_3d(0, b1, b2, neighbours_table1, neighbours_table2);
    time_step++;

    // End timer
    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};

    #pragma omp critical(io)
    {
        std::cout << "One sweep took " << elapsed_seconds.count() << " s\n";
        std::cout << "Expected execution time (not accounting for writing) =  " << 
            elapsed_seconds.count() * parameters.total_sweeps << " s " << std::endl;
    }

    for(int j = 1; j < 512; j++){
        do_cluster_sweep();
        update_observables(j, cluster_size);
        update_rg_observables_3d(j, b1, b2, neighbours_table1, neighbours_table2);
        time_step++;
    }

    #pragma omp critical (HDF5)
    {
        write_observables_after_loop(0);
        write_rg_observables_after_loop(0);
    }

    // End timer
    const auto finish512{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds512{finish512 - start};

    #pragma omp critical(io)
    {
        std::cout << "One loop of 512 sweeps took " << elapsed_seconds512.count() <<"s" << "\n";
        std::cout << "Need to do " << num_loops << " loops" << std::endl;
    }

    for (int i = 1; i < num_loops; i++){
        for(int j = 0; j < 512; j++){
            do_cluster_sweep();
            update_observables(j, cluster_size);
            update_rg_observables_3d(j, b1, b2, neighbours_table1, neighbours_table2);
            time_step++;
        }

        #pragma omp critical (HDF5)
        {
            write_observables_after_loop(i << 9);
            write_rg_observables_after_loop(i << 9);
        }
    }

    dvec final_x_magnetisation;
    dvec final_y_magnetisation;
    dvec final_energy;
    dvec final_average_cluster_size;
    dvec final_rg_x_magnetisation1;
    dvec final_rg_y_magnetisation1;
    dvec final_rg_x_magnetisation2;
    dvec final_rg_y_magnetisation2;
    dvec final_rg_x_magnetisation3{};
    dvec final_rg_y_magnetisation3{};
    dvec final_rg_energy1;
    dvec final_rg_energy2;
    dvec final_rg_energy3{};
    final_x_magnetisation.reserve(modulo);
    final_y_magnetisation.reserve(modulo);
    final_energy.reserve(modulo);
    final_average_cluster_size.reserve(modulo);
    final_rg_x_magnetisation1.reserve(modulo);
    final_rg_y_magnetisation1.reserve(modulo);
    final_rg_x_magnetisation2.reserve(modulo);
    final_rg_y_magnetisation2.reserve(modulo);
    final_rg_energy1.reserve(modulo);
    final_rg_energy2.reserve(modulo);

    for (int i = 0; i < modulo; i++){
        do_cluster_sweep();
        final_energy.emplace_back(model->compute_total_energy());
        final_x_magnetisation.emplace_back(model->compute_spin_magnetic_term(0));
        final_y_magnetisation.emplace_back(model->compute_spin_magnetic_term(1));
        final_average_cluster_size.emplace_back(compute_average_cluster_size(cluster_size));

        // RG 1
        model->compute_reduced_lattice_3d(b1, parameters.N, parameters.L);
        final_rg_x_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(0));
        final_rg_y_magnetisation1.emplace_back(model->compute_rg_spin_magnetic_term(1));
        final_rg_energy1.emplace_back(model->compute_rg_energy(b1, parameters.N, neighbours_table1, parameters.dim));
        // RG 2
        model->compute_reduced_lattice_3d(b2, parameters.N, parameters.L);
        final_rg_x_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(0));
        final_rg_y_magnetisation2.emplace_back(model->compute_rg_spin_magnetic_term(1));
        final_rg_energy2.emplace_back(model->compute_rg_energy(b2, parameters.N, neighbours_table2, parameters.dim));
        time_step++;
    }

    #pragma omp critical(HDF5)
    {
        write_observables_final(num_loops << 9, final_energy, final_x_magnetisation,final_y_magnetisation, final_average_cluster_size);
        write_rg_observables_final(num_loops << 9, final_rg_energy1, final_rg_energy2,final_rg_energy3, final_rg_x_magnetisation1, final_rg_x_magnetisation2,final_rg_x_magnetisation3, final_rg_y_magnetisation1, final_rg_y_magnetisation2, final_rg_y_magnetisation3);
        file->flush();
        close_writing();
    }

    #pragma omp critical(io)
    {
        std::cout << "Finished run for " << parameters.project_folder_path.filename() << "\n";
    }
}

void Simulation::do_metropolis_step() {
    double energy_diff = model->compute_energy_diff_flip();
    double r = rng::random_real_number(rng::engine);
    double probability = exp(-parameters.beta * energy_diff);
    if (r > probability) {
        model->revert_state_of_system();
    }
}

/** 
* Here time is measured as MCS/site, Monte Carlo Steps per site.
* So 1 step is N flips for N particles
*/
void Simulation::do_metropolis_sweep(){
    for (int i = 0; i < parameters.N; i++){
        do_metropolis_step();
    }
}

void Simulation::do_metropolis_recording_sweep() {
    // for (int i = 0; i < parameters.N; i++){
    //     do_metropolis_step();
    //     update_observables();
    // }
}

/**
* Here we need to normalise the MCS steps things (or not)
* I can do this like this at the moment and leave it.
*
* Total MCS = {\sum (size of each cluster)} / number particles
*/
void Simulation::do_cluster_step(){
    // Get the random stuff
    int random_index = rng::random_int_number(rng::engine, 0, parameters.N - 1);
    double angle_r = rng::random_angle(rng::engine) / 2;
    // Cannot do this yet
    model->flip_spin(random_index, angle_r);
    
    // Make stack to store indices used
    cluster_stack.clear();
    cluster_stack.emplace_back(random_index);
    std::fill(visited.begin(), visited.end(), 0);
    // cluster_stack.reserve(parameters.N);

    spins_flipped++;

    while(!cluster_stack.empty()){
        int current_index = cluster_stack.back();
        cluster_stack.pop_back();
        model->cluster_flip_neighbours(
            current_index,
            angle_r,
            cluster_stack,
            spins_flipped, 
            visited,
            parameters.dim
        );
        
    } 
}

void Simulation::do_cluster_sweep(){
    cluster_size.clear();
    spins_flipped = 0;
    int flip = 0;
        while (spins_flipped < parameters.N) {
            flip = spins_flipped;
            do_cluster_step();
            cluster_size.emplace_back(spins_flipped - flip);
        }
}

ivec Simulation::extract_int_vector_toml(const toml::array & array){
    ivec int_vector;
    int_vector.reserve(array.size());

    for (const auto &i : array) {
        int value = i.as_integer()->get();
        int_vector.emplace_back(value);
    }
    return int_vector;
}

dvec Simulation::extract_double_vector_toml(const toml::array& array) {
    dvec double_vector;
    double_vector.reserve(array.size());
    
    for (const auto &i: array) {
        double value = i.as_floating_point()->get();
        double_vector.emplace_back(value);
    }

    return double_vector;
}

void Simulation::update_observables(int position, ivec &cluster_size) {
    observables.energy_array[position] = model->compute_total_energy();
    observables.average_cluster_size[position] = compute_average_cluster_size(cluster_size);
    observables.x_magnetisation[position] = model->compute_spin_magnetic_term(0);
    observables.y_magnetisation[position] = model->compute_spin_magnetic_term(1);
}

void Simulation::update_rg_observables(int position, int b1, int b2,int b3, ivec & neigh_table1, ivec & neigh_table2, ivec & neigh_table3) {
    model->compute_reduced_lattice(b1, parameters.N, parameters.L);
    observables.rg_energy1[position] = model->compute_rg_energy(b1, parameters.N, neigh_table1, parameters.dim);
    observables.rg_x_magnetisation1[position] = model->compute_rg_spin_magnetic_term(0);
    observables.rg_y_magnetisation1[position] = model->compute_rg_spin_magnetic_term(1);

    model->compute_reduced_lattice(b2, parameters.N, parameters.L);
    observables.rg_energy2[position] = model->compute_rg_energy(b2, parameters.N, neigh_table2, parameters.dim);
    observables.rg_x_magnetisation2[position] = model->compute_rg_spin_magnetic_term(0);
    observables.rg_y_magnetisation2[position] = model->compute_rg_spin_magnetic_term(1);

    model->compute_reduced_lattice(b3, parameters.N, parameters.L);
    observables.rg_energy3[position] = model->compute_rg_energy(b3, parameters.N, neigh_table3, parameters.dim);
    observables.rg_x_magnetisation3[position] = model->compute_rg_spin_magnetic_term(0);
    observables.rg_y_magnetisation3[position] = model->compute_rg_spin_magnetic_term(1);
}

void Simulation::update_rg_observables_3d(int position, int b1, int b2, ivec & neigh_table1, ivec & neigh_table2) {
    model->compute_reduced_lattice_3d(b1, parameters.N, parameters.L);
    observables.rg_energy1[position] = model->compute_rg_energy(b1, parameters.N, neigh_table1, parameters.dim);
    observables.rg_x_magnetisation1[position] = model->compute_rg_spin_magnetic_term(0);
    observables.rg_y_magnetisation1[position] = model->compute_rg_spin_magnetic_term(1);

    model->compute_reduced_lattice_3d(b2, parameters.N, parameters.L);
    observables.rg_energy2[position] = model->compute_rg_energy(b2, parameters.N, neigh_table2, parameters.dim);
    observables.rg_x_magnetisation2[position] = model->compute_rg_spin_magnetic_term(0);
    observables.rg_y_magnetisation2[position] = model->compute_rg_spin_magnetic_term(1);
}

double Simulation::compute_average_cluster_size(ivec &cluster_size) {
    double average = 0;
    for (auto &i:cluster_size) {
       average+=i; 
    }

    return average / cluster_size.size();
}

void Simulation::write_lattice(int time) {
    model->write_lattice(parameters.lattice_set.get(),time, parameters.N);
}

void Simulation::write_observables_after_loop(hsize_t start) {
    assert(observables.energy_array.size() == 512 && "Energy array needs to be 512");
    assert(observables.average_cluster_size.size() == 512 && "Energy array needs to be 512");
    assert(observables.x_magnetisation.size() == 512 && "X_Magnetisation array needs to be 512");
    assert(observables.y_magnetisation.size() == 512 && "Y_Magnetisation array needs to be 512");
    parameters.magnetisation_set->select({start, 0}, {512, 1}).write_raw(observables.x_magnetisation.data(), HighFive::AtomicType<double>{});
    parameters.magnetisation_set->select({start, 1}, {512, 1}).write_raw(observables.y_magnetisation.data(), HighFive::AtomicType<double>{});;

    parameters.energy_set->select({start}, {512}).write(observables.energy_array);
    parameters.cluster_size_set->select({start}, {512}).write(observables.average_cluster_size);
}

void Simulation::write_rg_observables_after_loop(hsize_t start) {
    assert(observables.rg_energy1.size() == 512 && "Energy array needs to be 512");
    assert(observables.rg_energy2.size() == 512 && "Energy array needs to be 512");
    assert(observables.rg_x_magnetisation1.size() == 512 && "X_Magnetisation array needs to be 512");
    assert(observables.rg_x_magnetisation2.size() == 512 && "X_Magnetisation array needs to be 512");
    assert(observables.rg_y_magnetisation1.size() == 512 && "Y_Magnetisation array needs to be 512");
    assert(observables.rg_y_magnetisation2.size() == 512 && "Y_Magnetisation array needs to be 512");
    parameters.rg_magnetisation_set1->select({start, 0}, {512, 1}).write_raw(observables.rg_x_magnetisation1.data(), HighFive::AtomicType<double>{});
    parameters.rg_magnetisation_set1->select({start, 1}, {512, 1}).write_raw(observables.rg_y_magnetisation1.data(), HighFive::AtomicType<double>{});;
    parameters.rg_magnetisation_set2->select({start, 0}, {512, 1}).write_raw(observables.rg_x_magnetisation2.data(), HighFive::AtomicType<double>{});
    parameters.rg_magnetisation_set2->select({start, 1}, {512, 1}).write_raw(observables.rg_y_magnetisation2.data(), HighFive::AtomicType<double>{});;

    parameters.rg_energy_set1->select({start}, {512}).write(observables.rg_energy1);
    parameters.rg_energy_set2->select({start}, {512}).write(observables.rg_energy2);

    if (parameters.dim == 2){
        assert(observables.rg_energy3.size() == 512 && "Energy array needs to be 512");
        assert(observables.rg_x_magnetisation3.size() == 512 && "X_Magnetisation array needs to be 512");
        assert(observables.rg_y_magnetisation3.size() == 512 && "Y_Magnetisation array needs to be 512");
        parameters.rg_magnetisation_set3->select({start, 0}, {512, 1}).write_raw(observables.rg_x_magnetisation3.data(), HighFive::AtomicType<double>{});
        parameters.rg_magnetisation_set3->select({start, 1}, {512, 1}).write_raw(observables.rg_y_magnetisation3.data(), HighFive::AtomicType<double>{});;
        parameters.rg_energy_set3->select({start}, {512}).write(observables.rg_energy3);

    }
}

void Simulation::write_observables_final(hsize_t start, dvec energy, dvec x_magnetisation,dvec y_magnetisation, dvec cluster_size) {
    parameters.magnetisation_set->select({start, 0}, {x_magnetisation.size(), 1}).write_raw(x_magnetisation.data(), HighFive::AtomicType<double>{});
    parameters.magnetisation_set->select({start, 1}, {y_magnetisation.size(), 1}).write_raw(y_magnetisation.data(), HighFive::AtomicType<double>{});

    parameters.energy_set->select({start}, {energy.size()}).write(energy);
    parameters.cluster_size_set->select({start}, {cluster_size.size()}).write(cluster_size);
}
void Simulation::write_rg_observables_final(hsize_t start, dvec rg_energy1, dvec rg_energy2, dvec rg_energy3, dvec rg_x_magnetisation1,dvec rg_x_magnetisation2,dvec rg_x_magnetisation3, dvec rg_y_magnetisation1,dvec rg_y_magnetisation2, dvec rg_y_magnetisation3) {
    parameters.rg_magnetisation_set1->select({start, 0}, {rg_x_magnetisation1.size(), 1}).write_raw(rg_x_magnetisation1.data(), HighFive::AtomicType<double>{});
    parameters.rg_magnetisation_set1->select({start, 1}, {rg_y_magnetisation1.size(), 1}).write_raw(rg_y_magnetisation1.data(), HighFive::AtomicType<double>{});

    parameters.rg_magnetisation_set2->select({start, 0}, {rg_x_magnetisation2.size(), 1}).write_raw(rg_x_magnetisation2.data(), HighFive::AtomicType<double>{});
    parameters.rg_magnetisation_set2->select({start, 1}, {rg_y_magnetisation2.size(), 1}).write_raw(rg_y_magnetisation2.data(), HighFive::AtomicType<double>{});
    parameters.rg_energy_set1->select({start}, {rg_energy1.size()}).write(rg_energy1);
    parameters.rg_energy_set2->select({start}, {rg_energy2.size()}).write(rg_energy2);

    if (parameters.dim == 2){
        parameters.rg_magnetisation_set3->select({start, 0}, {rg_x_magnetisation3.size(), 1}).write_raw(rg_x_magnetisation3.data(), HighFive::AtomicType<double>{});
        parameters.rg_magnetisation_set3->select({start, 1}, {rg_y_magnetisation3.size(), 1}).write_raw(rg_y_magnetisation3.data(), HighFive::AtomicType<double>{});
    parameters.rg_energy_set3->select({start}, {rg_energy3.size()}).write(rg_energy3);
    }

}
