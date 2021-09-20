#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <DeepPhysicsSofa/Solver/LDLTSolver.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace DeepPhysicsSofa::solver::python {

    template <typename EigenMatrix>
    void bind_LDLTSolver(pybind11::module & m) {
        namespace py = pybind11;
        py::class_<LDLTSolver<EigenMatrix>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<LDLTSolver<EigenMatrix>>> c(m, "LDLTSolver");

        c.def("A", &LDLTSolver<EigenMatrix>::A);

        c.def("eigen_values", []()
        {
            using DenseMatrix = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic>;
            Eigen::SelfAdjointEigenSolver<DenseMatrix> eigen_solver;
            eigen_solver.compute(DenseMatrix(LDLTSolver<EigenMatrix>::A))
            return eigen_solver.eigenvalues().eval();
        })

        c.def("eigen_vectors", []()
        {
            using DenseMatrix = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic>;
            Eigen::SelfAdjointEigenSolver<DenseMatrix> eigen_solver;
            eigen_solver.compute(DenseMatrix(LDLTSolver<EigenMatrix>::A))
            return eigen_solver.eigenvectors().eval();
        })

        c.def("assemble", [](LDLTSolver<EigenMatrix> & solver, double m, double b, double k) {
            sofa::core::MechanicalParams mparams;
            mparams.setMFactor(m);
            mparams.setBFactor(b);
            mparams.setKFactor(k);
            solver.assemble(&mparams);
        }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));

        sofapython3::PythonFactory::registerType<LDLTSolver<EigenMatrix>>([](sofa::core::objectmodel::Base* o) {
            return py::cast(dynamic_cast<LDLTSolver<EigenMatrix>*>(o));
        });
    }

    void addLDLTSolver(pybind11::module &m);

}
