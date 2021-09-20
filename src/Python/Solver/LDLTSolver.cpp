#include "LDLTSolver.h"

#include <DeepPhysicsSofa/Solver/LDLTSolver.h>

namespace py = pybind11;

namespace SofaCaribou::solver::python {

    void addLDLTSolver(py::module & m) {
        bind_LDLTSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>(m);
    }

} // namespace SofaCaribou::solver::python