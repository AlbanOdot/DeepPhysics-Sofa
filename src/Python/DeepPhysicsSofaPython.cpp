#include <pybind11/pybind11.h>

#include <src/Python/Solver/LDLTSolver.h>

#include <vector>
#include <pybind11/stl_bind.h>

PYBIND11_MODULE(SofaCaribou, m) {
m.doc() = "DeepPhysicsSofa module";

// Solver bindings
DeepPhysicsSofa::solver::python::addLDLTSolver(m);

}