#pragma once

#include <SofaCaribou/Ode/StaticODESolver.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalVOpVisitor.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaCaribou/Solver/LinearSolver.h>
namespace DeepPhysicsSofa {

namespace ode {

    class HybridNewtonRaphson : public SofaCaribou::ode::StaticODESolver {
    public:
        SOFA_CLASS(HybridNewtonRaphson, StaticODESolver);
        HybridNewtonRaphson();

        template<typename T>
        using Data = sofa::core::objectmodel::Data<T>;

    protected:
        void solve(const sofa::core::ExecParams *params, SReal dt, sofa::core::MultiVecCoordId x_id,
                   sofa::core::MultiVecDerivId v_id) override;

        //Evaluate residual
        void computePredictionResidual(const sofa::core::ExecParams *params,
                                       sofa::core::MechanicalParams& mechanical_parameters,
                                       sofa::core::MultiVecDerivId& f_id,
                                       sofa::component::linearsolver::DefaultMultiMatrixAccessor& accessor,
                                       SofaCaribou::solver::LinearSolver * linear_solver);

        //Reset Simulation
        void clearSimulationData(sofa::simulation::common::VectorOperations &vop,
                             sofa::simulation::common::MechanicalOperations &mop,
                             sofa::core::MultiVecDerivId& f_id,
                             sofa::core::MultiVecDerivId& dx_id,
                             sofa::component::linearsolver::DefaultMultiMatrixAccessor& accessor,
                             SofaCaribou::solver::LinearSolver * linear_solver,
                             int newton_iterations);

        //Wheter or not we compute the tangeant stiffness matrix and inverse it
        Data<bool> d_solve_inversion;

        //Wheter or not we assemble the stiffness matrix before the PredictBeginEvent
        Data<bool> d_early_assemble;

        //Wheter or not we assemble the stiffness matrix before the PredictBeginEvent
        Data<bool> d_post_assemble;

        // Value of the residual from the prediction state
        Data<double> d_prediction_residual;

    };
} //ode
}// DeepPhysics
