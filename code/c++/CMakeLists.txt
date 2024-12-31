# Minimum required version of CMake
cmake_minimum_required(VERSION 3.14)

# Project name and language
project(SGEN LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the include directories
include_directories(include)

# Collect source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Create a static library
add_library(SGEN STATIC ${SOURCES})
set_target_properties(SGEN PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Add dependencies (optional)
find_package(Armadillo REQUIRED)
target_include_directories(SGEN PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(SGEN PUBLIC ${ARMADILLO_LIBRARIES})

# # Add an executable for the main .cpp file
# add_executable(optimised_moments main/optimised_moments.cpp)

# # Link the library to the executable
# target_link_libraries(optimised_moments PRIVATE SGEN)

# # Specify include directories for the executable
# target_include_directories(optimised_moments PRIVATE include)

# Add Python bindings
# Specify pybind11 directory if not found automatically
# set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "../extern/pybind11/build/")
find_package(pybind11 REQUIRED)

# Collect all binding files
set(BINDING_SOURCES
  bindings/bindings.cpp
  bindings/compartments/Compartment_binding.cpp
  bindings/compartments/Soma_binding.cpp
  bindings/compartments/Synapse_binding.cpp
  bindings/compartments/Dendritic_segment_binding.cpp
  # bindings/Neuron_binding.cpp
  bindings/engines/Analytic_engine_binding.cpp
)

# Add Python bindings module
pybind11_add_module(SGEN_Py ${BINDING_SOURCES})

# Link core library to the Python bindings
target_link_libraries(SGEN_Py PRIVATE SGEN)

# Include directories for Python bindings
target_include_directories(SGEN_Py PRIVATE ${CMAKE_SOURCE_DIR}/include)





















# # Minimum required version of CMake
# cmake_minimum_required(VERSION 3.14)

# # Project name and language
# project(SGEN LANGUAGES CXX)

# # Set the C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Specify the include directories
# include_directories(include)

# # Collect source files
# file(GLOB_RECURSE SOURCES "src/*.cpp")

# # Create a static library
# add_library(SGEN STATIC ${SOURCES})

# # Add dependencies (optional)
# find_package(Armadillo REQUIRED)
# target_include_directories(SGEN PUBLIC ${ARMADILLO_INCLUDE_DIRS})
# target_link_libraries(SGEN PUBLIC ${ARMADILLO_LIBRARIES})

# # Add an executable for the main .cpp file
# add_executable(optimised_moments main/optimised_moments.cpp)

# # Link the library to the executable
# target_link_libraries(optimised_moments PRIVATE SGEN)

# # Specify include directories for the executable
# target_include_directories(optimised_moments PRIVATE include)



# find_package(pybind11 REQUIRED)

# # Collect all binding files
# set(BINDING_SOURCES
#     bindings/bindings.cpp
#     bindings/neuron_bindings.cpp
#     bindings/compartments/soma_bindings.cpp
#     bindings/compartments/synapse_bindings.cpp
#     bindings/compartments/dendritic_segment_bindings.cpp
#     bindings/engines/analytic_engine_binding.cpp 
# )

# # Add Python bindings module
# pybind11_add_module(SGEN_Py ${BINDING_SOURCES})

# # Include directories for Python bindings
# target_include_directories(SGEN_Py PRIVATE ${CMAKE_SOURCE_DIR}/include)

# # Link core library
# target_link_libraries(SGEN_Py PRIVATE SGEN)





# # Add tests
# # option(BUILD_TESTS "Build tests" ON)
# # if (BUILD_TESTS)
# #     add_subdirectory(tests)
# # endif()

# # Installation instructions (optional)
# # install(TARGETS SGEN DESTINATION lib)
# # install(DIRECTORY include/ DESTINATION include)
