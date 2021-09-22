#include <src/Solver/LDLTSolver.h>
#include <SofaCaribou/Solver/LDLTSolver.inl>
DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

#ifdef DeepPhysicsSofa_WITH_MKL
    // Bug introduced in Eigen 3.3.8, fixed in bfdd4a9
    #ifndef EIGEN_USING_STD
        #define EIGEN_USING_STD(a) EIGEN_USING_STD_MATH(a)
    #endif
    #include <Eigen/PardisoSupport>
#endif


namespace DeepPhysicsSofa::solver {

    template<class EigenSolver_t>
    LDLTSolver<EigenSolver_t>::LDLTSolver()
            : d_compute_eigs(initData(&d_compute_eigs, false, "compute_eigs", "Recompute the eigen matrix and values during the solve"))
            , d_solve(initData(&d_solve, true, "solve", "True, actually solves the system, False : Only compute the factorisation and stuffs"))
    {}

    template<class EigenSolver_t>
    bool LDLTSolver<EigenSolver_t>::solve(const sofa::defaulttype::BaseVector * F,
                                          sofa::defaulttype::BaseVector *X) {

        bool solve_success = true;
        if(d_solve.getValue()) {
            solve_success = Base::solve(F, X);
        }

        if(d_compute_eigs.getValue())
        {
            p_eigen_solver.compute(this->A()->matrix());
            p_eigen_values = p_eigen_solver.eigenvalues().eval();
            p_eigen_vectors = p_eigen_solver.eigenvectors().eval();
        }
        return solve_success;
    }


    static int DPSSparseLDLTSolverClass = sofa::core::RegisterObject("DeepPhysicsSofa Sparse LDLT linear solver")
                                               .add< LDLTSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>> >(true)
#ifdef DeepPhysicsSofa_WITH_MKL
    .add< LDLTSolver<Eigen::PardisoLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>> >()
#endif
    ;

} // namespace DeepPhysicsSofa::solver