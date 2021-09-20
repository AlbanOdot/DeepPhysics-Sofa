project(DeepPhysicsSofa.Python)

set(HEADER_FILES
        Solver/LDLTSolver.h
        )

set(SOURCE_FILES
        DPS.cpp
        Solver/LDLTSolver.cpp
        )

set(PYTHON_FILES
        __init__.py
        )


find_package(SofaPython3 REQUIRED)
DPS_add_python_module(DeepPhysicsSofa
        TARGET_NAME       ${PROJECT_NAME}
        TARGET_ALIAS      DeepPhysicsSofa::Python.DeepPhysicsSofa
        TARGET_DEPENDS    Caribou::Geometry Caribou::Python.Config DeepPhysicsSofa SofaCaribou SofaPython3::Bindings.Sofa.Core
        COMPONENT_NAME    DeepPhysicsSofa
        PACKAGE_NAME      DeepPhysicsSofa
        DESTINATION       DeepPhysicsSofa
        SOURCE_FILES      ${SOURCE_FILES}
        PUBLIC_HEADERS    ${HEADER_FILES}
        PYTHON_FILES      ${PYTHON_FILES}
        PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
        )

list(APPEND target_rpath
        "$ORIGIN/../../../../../../lib"
        "$loader_path/../../../../../../lib"
        "$ORIGIN/../../../../lib"
        "$loader_path/../../../../lib"
        )
set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        INSTALL_RPATH "${target_rpath}"
)

if (SOFA_VERSION VERSION_GREATER_EQUAL "20.12.99")
    # Set RPATH towards relocatable dependencies
    sofa_auto_set_target_rpath(TARGETS ${PROJECT_NAME} RELOCATABLE "plugins")
endif()