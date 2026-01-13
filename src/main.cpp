#include "simulation/simulation.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <memory>
#include <toml++/toml.hpp>

/*****************************************
 * What is the algorithm process?
 *  1. get from argv location of toml file
 *  2. get parameters from toml
 *      - create a struct where to store all of them
 *  3. create h5 file
 *  4. select model with toml file
 *  5. initialise matrices and models
 *  6. select algorithm
 *  7. run algorithm
 *  8. write while the algorithm is running 
 *      beware with memory!!!
 *  9. close h5 file and exit simulation
 *  ! Error handling, assertions
 ****************************************/
int main(int argc, char *argv[]){
    std::string project_folder;
    std::string running_model;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg == "-p" && i + 1 < argc) {
                project_folder = argv[i + 1];
                ++i; // Skip next argument since it's the folder name
            }
            else if (arg == "-m" && i + 1 < argc){
                running_model = argv[i + 1];
                ++i;
            }
        }

        if (!project_folder.empty()) {
            // Assuming it's being executed from build folder
            std::filesystem::path project_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
            std::cout << "Using project folder: " << project_path << "\n";
        }     

    
    try{ 
        // Iterate over the different directories
        std::filesystem::path project_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
        std::vector<std::filesystem::path> directories;
        // Separate this for easier parallelisation implementation
        for (auto const& dir_entry : std::filesystem::directory_iterator{project_path}){ 
            std::cout <<"dir entry path"<< dir_entry.path() << '\n';
            directories.push_back(dir_entry.path());
        }
        for (int i = 0; i < directories.size(); i++){
            // Using smart pointers
            std::unique_ptr<Simulation> sim = std::make_unique<Simulation>();
            sim->parse_parameters(directories[i], running_model);
            sim->run();
        }
        return 0;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Configuration error: " << e.what() << std::endl;
        return 4;
    }
}
