#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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
            std::cout << "Using project folder: ./" << projectFolder << std::endl;
        }     


}
