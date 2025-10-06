cmake_minimum_required(VERSION 3.13)
project(capstone)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/build
)
add_subdirectory(src)
