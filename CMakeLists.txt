cmake_minimum_required(VERSION 3.12)
project(DeepPhysicsSofa VERSION 1.0.0)

# Policies
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0094 NEW)

set(DeepPhysicsSofa_VERSION ${PROJECT_VERSION})

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

# Detect if DeepPhysicsSofa is a subproject of SOFA
if ("${CMAKE_PROJECT_NAME}" STREQUAL "Sofa")
    set(DeepPhysicsSofa_COMPILED_AS_SOFA_SUBPROJECT 1)
    message(STATUS "${PROJECT_NAME} is compiled as a subproject to SOFA.")
endif ()

if (DeepPhysicsSofa_COMPILED_AS_SOFA_SUBPROJECT)
    option(PLUGIN_DeepPhysicsSofa "Build the DeepPhysicsSofa plugin." ON)
    if (NOT PLUGIN_DeepPhysicsSofa)
        return()
    else()
        message("Adding plugin DeepPhysicsSofa")
    endif()
endif()

# Include CMake specific tools
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake")
list(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_BINARY_DIR}/cmake")
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(DeepPhysicsSofaTools)
set(CMAKE_REQUIRED_QUIET 1)


CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

# OPTIONS
option(DeepPhysicsSofa_USE_FLOAT "Use float precision for floating point values." OFF)
CMAKE_DEPENDENT_OPTION(DeepPhysicsSofa_OPTIMIZE_FOR_NATIVE "Build with -march=native to allow msse support" ON "COMPILER_SUPPORTS_MARCH_NATIVE" OFF)

# SOFA option
find_package(SOFA COMPONENTS SofaFramework QUIET)
CMAKE_DEPENDENT_OPTION(DeepPhysicsSofa_WITH_SOFA "Compile the SofaDeepPhysicsSofa library" ON "SOFA_FOUND" OFF)
if (DeepPhysicsSofa_WITH_SOFA)
    list(APPEND CMAKE_PREFIX_PATH ${SOFA_MODULE_DIR})
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    find_package(SofaPython3 QUIET) # We need to do this first to setup the exact python version required
endif()

if (DeepPhysicsSofa_WITH_SOFA)
    get_filename_component(SOFA_ROOT_DIR "${SofaFramework_DIR}/../../.." ABSOLUTE)
    message(STATUS
            "DeepPhysicsSofa with SOFA support\n"
            "    SOFA version: ${SofaFramework_VERSION}\n"
            "    SOFA location: ${SOFA_ROOT_DIR}"
            )
endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(NOT WIN32)
    if ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*[cC][lL][aA][nN][gG]")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

add_subdirectory(src)

install(
        FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeepPhysicsSofaTools.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSOFA.cmake"
        DESTINATION lib/cmake/DeepPhysicsSofa
        COMPONENT headers
)
