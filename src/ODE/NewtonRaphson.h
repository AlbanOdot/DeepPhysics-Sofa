#pragma once

#include <SofaCaribou/config.h>
#include <src/ODE/ODEInterface.h>

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
 * This class implements a NewtonRaphson solver for SOFA.
 */
    class NewtonRaphson : public ODEInterface {
    public:
        SOFA_CLASS(NewtonRaphson, ODEInterface);

        template <typename T>
        using Data = sofa::core::objectmodel::Data<T>;

        template <typename T>
        using Link = sofa::core::objectmodel::SingleLink<NewtonRaphson, T, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

        CARIBOU_API
        NewtonRaphson();

        CARIBOU_API
        void solve (const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id) override;

    protected:

        virtual void assemble_rhs_vector(const sofa::core::MechanicalParams & mechanical_parameters,
                                         const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                         sofa::core::MultiVecDerivId & f_id,
                                         sofa::defaulttype::BaseVector * f);

        virtual void assemble_system_matrix(const sofa::core::MechanicalParams & mechanical_parameters,
                                            sofa::component::linearsolver::DefaultMultiMatrixAccessor & matrix_accessor,
                                            sofa::defaulttype::BaseMatrix * A);

        virtual void propagate_solution_increment(const sofa::core::MechanicalParams & mechanical_parameters,
                                                  const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                                  const sofa::defaulttype::BaseVector * dx,
                                                  sofa::core::MultiVecCoordId & x_id,
                                                  sofa::core::MultiVecDerivId & v_id,
                                                  sofa::core::MultiVecDerivId & dx_id);
    };
}
