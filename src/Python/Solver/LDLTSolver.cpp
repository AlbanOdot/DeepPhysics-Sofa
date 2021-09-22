#include "LDLTSolver.h"

#include <src/Solver/LDLTSolver.h>
//#ifdef DeepPhysicsSofa_WITH_MKL
//// Bug introduced in Eigen 3.3.8, fixed in bfdd4a9
//    #ifndef EIGEN_USING_STD
//        #define EIGEN_USING_STD(a) EIGEN_USING_STD_MATH(a)
//    #endif
//    #include <Eigen/PardisoSupport>
//#endif
namespace py = pybind11;

namespace DeepPhysicsSofa::solver::python {

    void addLDLTSolver(py::module & m) {
        bind_LDLTSolver< Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>,Eigen::Lower, Eigen::AMDOrdering<int>>>(m);
//        bind_LDLTSolver<Eigen::PardisoLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>>(m);
    }

} // namespace DeepPhysicsSofa::solver::python