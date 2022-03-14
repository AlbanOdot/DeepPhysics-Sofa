//
// Created by alban on 14/03/2022.
//

#include <src/ODE/KSolver.h>
#include <iomanip>
#include <chrono>

#include <SofaCaribou/Solver/LinearSolver.h>
#include <SofaCaribou/Solver/LDLTSolver.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>
#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

#include <src/Event/PredictBeginEvent.h>
#include <src/Event/PredictEndEvent.h>
#include <src/Event/PredictionPickedEvent.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/PropagateEventVisitor.h>

namespace DeepPhysicsSofa {
    namespace ode {
        KSolver::KSolver() {}

        void KSolver::solve(const sofa::core::ExecParams *params, SReal dt, sofa::core::MultiVecCoordId x_id,
                            sofa::core::MultiVecDerivId v_id)
        {
            using namespace sofa::helper::logging;
            using namespace std::chrono;
            using std::chrono::steady_clock;

            // Make sure we have a linear solver, and that it implements the Caribou::solver::LinearSolver interface
            static bool error_message_already_printed = false;
            if (not has_valid_linear_solver()) {
                if (not error_message_already_printed) {
                    msg_error() << "The system will NOT be solved.";
                    msg_error() << "No compatible linear solver have been set. Use the '"
                                << l_linear_solver.getName()
                                << "' attribute to specify the path towards a linear solver.";
                    error_message_already_printed = true;
                }
                return;
            }
            error_message_already_printed = false;

            // Get the current context
            const auto context = this->getContext();

            // Set the multi-vector identifier inside the mechanical parameters.
            sofa::core::MechanicalParams mechanical_parameters (*params);
            mechanical_parameters.setX(x_id);
            mechanical_parameters.setF(sofa::core::ConstVecDerivId::force());
            mechanical_parameters.setDf(sofa::core::ConstVecDerivId::dforce());
            mechanical_parameters.setDt(dt);

            // Get the linear solver that implements the SofaCaribou::solver::LinearSolver interface
            auto linear_solver = dynamic_cast<SofaCaribou::solver::LinearSolver *>(l_linear_solver.get());

            // Create the vector and mechanical operations tools. These are used to execute special operations (multiplication,
            // additions, etc.) on multi-vectors (a vector that is stored in different buffers inside the mechanical objects)
            sofa::simulation::common::VectorOperations vop( &mechanical_parameters, context );
            sofa::simulation::common::MechanicalOperations mop( &mechanical_parameters, context );

            // Let the mechanical operations know that this is an implicit solver. This will be propagated back to the
            // force fields during the addForce and addKToMatrix phase, which will let them recompute their internal
            // stresses if they have a non-linear relationship with the displacement.
            mop->setImplicit(true);

            // Right hand side term (internal + external forces)
            auto f_id = sofa::core::MultiVecDerivId(sofa::core::VecDerivId::force());
            vop.v_clear(f_id);

            // Start the advanced timer
            sofa::helper::ScopedAdvancedTimer timer (this->getClassName() + "::Solve");
            // ###########################################################################
            // #                           Mechanical graph                              #
            // ###########################################################################
            // # Construct the mechanical graph by finding top level mechanical objects, #
            // # mechanical mapping, and mapped mechanical objects. This graph will be   #
            // # used to compute the final assembled system matrix.                      #
            // ###########################################################################

            // For now, let the "default" multi-matrix accessor go down the scene graph and
            // accumulate the mechanical objects and mappings. This one will not really
            // compute the mechanical graph (not explicitly at least). Hence the following
            sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

            // Step 1   Get dimension of each top level mechanical states using
            //          BaseMechanicalState::getMatrixSize(), and accumulate mechanical
            //          objects and mapping matrices
            mop.getMatrixDimension(nullptr, nullptr, &accessor);
            const auto n = static_cast<sofa::Size>(accessor.getGlobalDimension());

            // Step 2   Does nothing more than to accumulate from the previous step a list of
            //          "MatrixRef = <MechanicalState*, MatrixIndex>" where MatrixIndex is the
            //          (i,i) position of the given top level MechanicalState* inside the global
            //          system matrix. This global matrix hence contains one sub-matrix per top
            //          level mechanical state.
            accessor.setupMatrices();
            // Step 3   Let the linear solver create the system matrix and vector buffers
            //          using the previously computed system size n
            p_A.reset(linear_solver->create_new_matrix(n, n));
            p_A->clear();
            p_F.reset(linear_solver->create_new_vector(n));
            p_F->clear();


            // ###########################################################################
            // #                             First residual                              #
            // ###########################################################################
            // # Before starting any newton iterations, we first need to compute         #
            // # the residual with the updated right-hand side (the new load increment)  #
            // ###########################################################################
            // If working with rest shapes, initial residual only contains external forces
            // We now call predictBeginEvent in order to catch the external forces
            PredictBeginEvent PBev ( this->getContext()->getRootContext()->getDt ());
            sofa::simulation::PropagateEventVisitor pTBev ( params, &PBev );
            this->getContext()->getRootContext()->executeVisitor(&pTBev );

            // Step 1   Assemble the force vector
            sofa::helper::AdvancedTimer::stepBegin("ComputeForce");
            this->assemble_rhs_vector(mechanical_parameters, accessor, f_id, p_F.get());
            sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

            // Step 2   Compute the initial residual
            d_prediction_residual.setValue(SofaCaribou::Algebra::dot(p_F.get(), p_F.get()));

            sofa::helper::AdvancedTimer::stepBegin("MBKBuild");
            p_A->clear();
            this->assemble_system_matrix(mechanical_parameters, accessor, p_A.get());
            linear_solver->set_system_matrix(p_A.get());
            sofa::helper::AdvancedTimer::stepEnd("MBKBuild");

            // If working with rest shapes, initial residual only contains external forces
            // We now call predictBeginEvent in order to catch the external forces
            PredictEndEvent PEev ( this->getContext()->getRootContext()->getDt ());
            sofa::simulation::PropagateEventVisitor pTPEev ( params, &PEev );
            this->getContext()->getRootContext()->executeVisitor(&pTPEev );
        }

        int KSolverClass = sofa::core::RegisterObject("K Solver").add< KSolver >();
    }
}