
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(Simulation)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
# CRY Files
#include_directories(/home/subhendu/Tools/cry_v1.7/src/)
#link_directories(/home/subhendu/Tools/cry_v1.7/lib/)

add_executable(sim sim.cc ${sources})#here add the sim or main file name

# target_link_libraries(main ${Geant4_LIBRARIES} CRY)
target_link_libraries(sim ${Geant4_LIBRARIES})

add_custom_target(Simulation DEPENDS sim)

set(SIM_SCRIPTS
    init_vis.mac
    vis.mac
    run.mac
)

foreach(_script ${SIM_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
    )
endforeach()
