# Minimum required version of CMake
cmake_minimum_required(VERSION 3.14)

# Project name and language
project(NeuronLib LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the include directories
include_directories(include)

# Collect source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Create a static library
add_library(neuron_lib STATIC ${SOURCES})

# Add dependencies (optional)
find_package(Armadillo REQUIRED)
target_include_directories(neuron_lib PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(neuron_lib PUBLIC ${ARMADILLO_LIBRARIES})

# Add an executable for the main .cpp file
add_executable(optimised_moments main/optimised_moments.cpp)

# Link the library to the executable
target_link_libraries(optimised_moments PRIVATE neuron_lib)

# Specify include directories for the executable
target_include_directories(optimised_moments PRIVATE include)

# Add tests
# option(BUILD_TESTS "Build tests" ON)
# if (BUILD_TESTS)
#     add_subdirectory(tests)
# endif()

# Installation instructions (optional)
# install(TARGETS neuron_lib DESTINATION lib)
# install(DIRECTORY include/ DESTINATION include)
