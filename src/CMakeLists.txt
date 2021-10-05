set(HEADER_FILES
        Event/PredictBeginEvent.h
        Event/PredictEndEvent.h
        Event/PredictionPickedEvent.h
        Forcefield/FreeConstantForceField.h
        ODE/HybridNewtonRaphson.h
        )

set(SOURCE_FILES
        Event/PredictBeginEvent.cpp
        Event/PredictEndEvent.cpp
        Event/PredictionPickedEvent.cpp
        Forcefield/FreeConstantForceField.cpp
        ODE/HybridNewtonRaphson.cpp
        )

# Caribou
find_package(Caribou COMPONENTS Algebra REQUIRED QUIET)

# Sofa's packages
find_package(SOFA COMPONENTS SofaFramework SofaBoundaryCondition SofaBaseLinearSolver SofaBaseTopology SofaEigen2Solver QUIET REQUIRED)
find_package(SofaCaribou)
# Bugged dependencies
if (macos)
    find_package(Boost COMPONENTS system QUIET) # Bug introduced in SOFA v21.06 where this dependencies isn't propagated
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})
add_library(DeepPhysicsSofa::Plugin ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBoundaryCondition SofaBaseTopology SofaBaseLinearSolver SofaEigen2Solver)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Algebra Caribou::Geometry Caribou::Topology Caribou::Mechanics SofaCaribou)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

list(APPEND target_rpath
        "$ORIGIN/../lib"
        "@executable_path/../lib"
        )

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")


DeepPhysicsSofa_add_target_to_component (
        TARGET_NAME        ${PROJECT_NAME}
        COMPONENT_NAME     DeepPhysicsSofa
        PACKAGE_NAME       DeepPhysicsSofa
        PUBLIC_HEADERS     "${TEMPLATE_FILES};${HEADER_FILES}"
        HEADER_SRC_PREFIX  "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

install (
        EXPORT DeepPhysicsSofaTargets
        DESTINATION "lib/cmake/DeepPhysicsSofa"
        COMPONENT headers
)

# This is required to be able to do "find_package(DeepPhysicsSofa)"
DeepPhysicsSofa_create_package(DeepPhysicsSofa ${DeepPhysicsSofa_VERSION})
# Make the target relocatable
set_target_properties(${PROJECT_NAME} PROPERTIES RELOCATABLE_INSTALL_DIR "plugins/DeepPhysicsSofa")
set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_PROPERTIES "RELOCATABLE_INSTALL_DIR")

if (SOFA_VERSION VERSION_GREATER_EQUAL "20.12.00")
    # Set RPATH towards relocatable dependencies
    sofa_auto_set_target_rpath(TARGETS ${PROJECT_NAME} RELOCATABLE "plugins")
endif()