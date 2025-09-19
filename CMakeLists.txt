cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
option(GEANT4_USE_QT "Enable Qt visualization" ON)

project(Simulation)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

# CRY includes and library
include_directories(/home/oem/cry_v1.7/src)
link_directories(/home/oem/cry_v1.7)

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

add_executable(sim sim.cc ${sources} ${headers})

# Link CRY as well
target_link_libraries(sim ${Geant4_LIBRARIES} /home/oem/MSC_DISSERTATION/15_08_25/Geometry_Build_4/cry_v1.7/lib/libCRY.a)

add_custom_target(Simulation DEPENDS sim)

