set(HEADER_FILES
        Event/PredictBeginEvent.h
        Event/PredictEndEvent.h
        ODE/HybridNewtonRaphson.h
        ODE/NewtonRaphson.h
        Solver/LDLTSolver.h)

set(SOURCE_FILES
        Event/PredictBeginEvent.cpp
        Event/PredictEndEvent.cpp
        ODE/HybridNewtonRaphson.cpp
        ODE/NewtonRaphson.cpp
        Solver/LDLTSolver.cpp
        )

# Caribou
find_package(Caribou COMPONENTS Algebra Geometry Topology Mechanics REQUIRED QUIET)

# Sofa's packages
find_package(SOFA COMPONENTS SofaFramework SofaBaseLinearSolver SofaBaseTopology SofaEigen2Solver QUIET REQUIRED)
find_package(SofaCaribou)
# Optional packages
find_package(Eigen3 QUIET REQUIRED)

# Bugged dependencies
if (macos)
    find_package(Boost COMPONENTS system QUIET) # Bug introduced in SOFA v21.06 where this dependencies isn't propagated
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})
add_library(DPS::Plugin ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBaseTopology SofaBaseLinearSolver SofaEigen2Solver)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Algebra Caribou::Geometry Caribou::Topology Caribou::Mechanics SofaCaribou)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

list(APPEND target_rpath
        "$ORIGIN/../lib"
        "@executable_path/../lib"
        )

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")


DPS_add_target_to_component (
        TARGET_NAME        ${PROJECT_NAME}
        COMPONENT_NAME     SofaCaribou
        PACKAGE_NAME       SofaCaribou
        PUBLIC_HEADERS     "${TEMPLATE_FILES};${HEADER_FILES}"
        HEADER_SRC_PREFIX  "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

find_package(SofaPython3 QUIET)
#add_subdirectory(Python)

#install (
#        EXPORT DPSTargets
#        DESTINATION "lib/cmake/DPS"
#        COMPONENT headers
#)

# This is required to be able to do "find_package(SofaCaribou)"
DPS_create_package(DPS ${DPS_VERSION})
# Make the target relocatable
set_target_properties(${PROJECT_NAME} PROPERTIES RELOCATABLE_INSTALL_DIR "plugins/DPS")
set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_PROPERTIES "RELOCATABLE_INSTALL_DIR")

if (SOFA_VERSION VERSION_GREATER_EQUAL "20.12.00")
    # Set RPATH towards relocatable dependencies
    sofa_auto_set_target_rpath(TARGETS ${PROJECT_NAME} RELOCATABLE "plugins")
endif()