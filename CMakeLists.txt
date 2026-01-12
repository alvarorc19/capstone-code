cmake_minimum_required(VERSION 3.13)
project(capstone)

## Custom targets
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/build
)


## Dependencies
include(FetchContent)

# find_package(HDF5 REQUIRED COMPONENTS CXX)
# link_directories(${HDF5_LIBRARY_DIRS})
# include_directories(${HDF5_INCLUDE_DIRS})

# find_package(Boost 1.74.0 REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

# find_package(OpenMP 4.5 REQUIRED COMPONENTS CXX)
# include_directories(${OpenMP_CXX_INCLUDE_DIRS})

FetchContent_Declare(
        highfive
        GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
#        GIT_REPOSITORY git@github.com:BlueBrain/HighFive.git
        GIT_SHALLOW ON
        GIT_TAG v3.2.0
)
FetchContent_GetProperties(highfive)
FetchContent_MakeAvailable(highfive)
include_directories(BEFORE SYSTEM ${highfive_SOURCE_DIR}/include ${highfive_BINARY_DIR}/include)


# FetchContent_Declare(
#         tomlplusplus
#         GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
#         GIT_TAG v3.4.0
# )
# FetchContent_GetProperties(tomlplusplus)
# FetchContent_MakeAvailable(tomlplusplus)
# include_directories(BEFORE SYSTEM ${tomlplusplus_SOURCE_DIR}/include ${tomlplusplus_BINARY_DIR}/include)

## Subdirectories
add_subdirectory(src)
