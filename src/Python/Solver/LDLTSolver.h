#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <src/Solver/LDLTSolver.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>


namespace DeepPhysicsSofa::solver::python {

    template <typename EigenSolver_t>
    void bind_LDLTSolver(pybind11::module & m) {
        namespace py = pybind11;
        using SOLVER = DeepPhysicsSofa::solver::LDLTSolver<EigenSolver_t>;
        py::class_<SOLVER, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<SOLVER>> c(m, "DPSLDLTSolver");

        c.def("A", &SOLVER::A);
        c.def("eigen_values", &SOLVER::eigen_values);
        c.def("eigen_vectors", &SOLVER::eigen_vectors);
        c.def("assemble", [](SOLVER & solver, double m, double b, double k) {
            sofa::core::MechanicalParams mparams;
            mparams.setMFactor(m);
            mparams.setBFactor(b);
            mparams.setKFactor(k);
            solver.assemble(&mparams);
        }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));
        sofapython3::PythonFactory::registerType<SOLVER>([](sofa::core::objectmodel::Base* o) {
            return py::cast(dynamic_cast<SOLVER*>(o));
        });
    }


    void addLDLTSolver(pybind11::module &m);

}
