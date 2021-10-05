#pragma once

#include <SofaCaribou/Ode/StaticODESolver.h>

namespace DeepPhysicsSofa::ode{

class HybridNewtonRaphson : public SofaCaribou::ode::StaticODESolver
{
public:
    SOFA_CLASS(HybridNewtonRaphson, StaticODESolver);
    HybridNewtonRaphson();

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

protected:
    void solve (const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id) override;
    // Only compute the residual of the system
    void solveResidual(const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id);

    // Compute both residual and ground truth of the system
    void solveInversion(const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id);

    //Wheter or not we compute the tangeant stiffness matrix and inverse it
    Data<bool> d_solve_inversion;

    // Value of the residual from the prediction state
    Data<double> d_prediction_residual;

    /// Global system RHS vector (the forces)
    std::unique_ptr<sofa::defaulttype::BaseVector> p_F_prediction;

    //Function pointer toward the corresponding way of solving. (usually "solveRasidual" during learning phase and "solveInversion" during testing/using phase)
    void (HybridNewtonRaphson::*m_solve)(const sofa::core::ExecParams* params /* PARAMS FIRST */, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult);

};

}