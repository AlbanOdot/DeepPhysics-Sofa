#include <pybind11/pybind11.h>

#include <DeepPhysicsSofa/Python/Solver/LDLTSolver.h>

#include <vector>
#include <pybind11/stl_bind.h>

PYBIND11_MODULE(SofaCaribou, m) {
m.doc() = "SofaCaribou module";

// Solver bindings
SofaCaribou::solver::python::addConjugateGradientSolver(m);

}