#pragma once

#include <SofaCaribou/config.h>
#include <src/ODE/NewtonRaphson.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/BaseVector.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/helper/OptionsGroup.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
        DISABLE_ALL_WARNINGS_END

#include <memory>

namespace DeepPhysicsSofa::ode {

/**
 * This class implements a Hybrid Newton Raphson solver for SOFA.
 *
 * Following this article https://hal.archives-ouvertes.fr/hal-03327818/document
 *
 */
    class HybridNewtonRaphson : public NewtonRaphson {
    public:
        SOFA_CLASS(HybridNewtonRaphson, NewtonRaphson);

        template <typename T>
        using Data = sofa::core::objectmodel::Data<T>;

        template <typename T>
        using Link = sofa::core::objectmodel::SingleLink<HybridNewtonRaphson, T, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

        CARIBOU_API
        HybridNewtonRaphson();

        CARIBOU_API
        void solve (const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id) override;

    protected:

        // Only compute the residual of the system
        void solveResidual(const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id);

        // Compute both residual and ground truth of the system
        void solveInversion(const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id);

        //Wheter or not we compute the tangeant stiffness matrix and inverse it
        Data<bool> d_solve_inversion;

        // Value of the residual from the prediction state
        Data<double> d_prediction_residual;

        //Function pointer toward the corresponding way of solving. (usually "solveRasidual" during learning phase and "solveInversion" during testing/using phase)
        void (HybridNewtonRaphson::*m_solve)(const sofa::core::ExecParams* params /* PARAMS FIRST */, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult);
    };
}
