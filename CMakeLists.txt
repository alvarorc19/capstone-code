cmake_minimum_required(VERSION 3.13)
project(capstone)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_EXTENSIONS OFF)

## Custom targets
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/build
)


## Dependencies
include(FetchContent)

find_package(HDF5 REQUIRED COMPONENTS CXX)
link_directories(${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

# find_package(Boost 1.74.0 REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

find_package(OpenMP REQUIRED COMPONENTS CXX)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()

include_directories(${OpenMP_CXX_INCLUDE_DIRS})

# FetchContent_Declare(
#         highfive
# 	#        GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
#         GIT_REPOSITORY git@github.com:BlueBrain/HighFive.git
#         GIT_SHALLOW ON
#         GIT_TAG v2.10.1
# )
# FetchContent_GetProperties(highfive)
# FetchContent_MakeAvailable(highfive)
find_package(highfive)
include_directories(BEFORE SYSTEM ${highfive_SOURCE_DIR}/include ${highfive_BINARY_DIR}/include)


FetchContent_Declare(
        tomlplusplus
        GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
        GIT_TAG v3.4.0
)
FetchContent_GetProperties(tomlplusplus)
FetchContent_MakeAvailable(tomlplusplus)
include_directories(BEFORE SYSTEM ${tomlplusplus_SOURCE_DIR}/include ${tomlplusplus_BINARY_DIR}/include)

## Subdirectories
add_subdirectory(src)
