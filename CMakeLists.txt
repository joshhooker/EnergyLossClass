#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ENERGYLOSS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(output test.cpp ${sources} ${headers})

# Move SRIM file to build directory
set(SCRIPTS
  CarbonCsI.dat
  CarbonCsIReduced.dat
  EnergyLoss.dat
  ProtonD2.dat
)

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
    )
endforeach()


#
#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ENERGYLOSS DEPENDS output)

