#pragma once

#include <SofaCaribou/Solver/LDLTSolver.h>
#include <Eigen/Eigenvalues>

namespace DeepPhysicsSofa::solver {

    template <class EigenSolver_t>
    class LDLTSolver : public SofaCaribou::solver::LDLTSolver< EigenSolver_t> {
        public:
            SOFA_CLASS(SOFA_TEMPLATE(LDLTSolver, EigenSolver_t), SOFA_TEMPLATE(SofaCaribou::solver::LDLTSolver, EigenSolver_t));

            template <typename T>
            using Data = sofa::Data<T>;

            using Base = SofaCaribou::solver::LDLTSolver< EigenSolver_t>;
            using EigenValuesSolver = Eigen::SelfAdjointEigenSolver<typename EigenSolver_t::MatrixType>;

            LDLTSolver();

            bool solve(const sofa::defaulttype::BaseVector * F, sofa::defaulttype::BaseVector * X) override;

            auto eigen_values () const -> const typename EigenValuesSolver::RealVectorType & {
                return p_eigen_values;
            }

            auto eigen_vectors ()  const -> const typename EigenValuesSolver::EigenvectorsType & {
                return p_eigen_vectors;
            }

            private:

                EigenValuesSolver p_eigen_solver;
                typename EigenValuesSolver::EigenvectorsType p_eigen_vectors;
                typename EigenValuesSolver::RealVectorType p_eigen_values;
                Data<bool> d_compute_eigs;
                Data<bool> d_solve;
            };
}