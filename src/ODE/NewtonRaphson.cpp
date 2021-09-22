#include <src/ODE/NewtonRaphson.h>

#include <iomanip>
#include <chrono>

#include <SofaCaribou/Solver/LinearSolver.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>
#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Size = int; }
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#else
#include <sofa/simulation/mechanicalvisitor/MechanicalApplyConstraintsVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalComputeForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorFromBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorToBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalResetForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalVOpVisitor.h>
using namespace sofa::simulation::mechanicalvisitor;
#endif
DISABLE_ALL_WARNINGS_END


namespace DeepPhysicsSofa::ode {

    using sofa::core::ExecParams;
    using sofa::core::MultiVecCoordId;
    using sofa::core::MultiVecDerivId;
    using namespace sofa::simulation;
    using Timer = sofa::helper::AdvancedTimer;

    NewtonRaphson::NewtonRaphson()
    {
        d_pattern_analysis_strategy.setValue(sofa::helper::OptionsGroup(std::vector < std::string > {
                "NEVER", "BEGINNING_OF_THE_SIMULATION", "BEGINNING_OF_THE_TIME_STEP", "ALWAYS"
        }));

        // Select the default value
        set_pattern_analysis_strategy(PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP);
    }

    void NewtonRaphson::solve(const ExecParams *params, SReal dt, MultiVecCoordId x_id, MultiVecDerivId v_id) {
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
        mechanical_parameters.setV(v_id);
        mechanical_parameters.setF(sofa::core::ConstVecDerivId::force());
        mechanical_parameters.setDf(sofa::core::ConstVecDerivId::dforce());
        mechanical_parameters.setDx(sofa::core::ConstVecDerivId::dx());
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

        // Options for the Newton-Raphson
        const auto   pattern_strategy = pattern_analysis_strategy();
        const auto & correction_tolerance_threshold = d_correction_tolerance_threshold.getValue();
        const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
        const auto & absolute_residual_tolerance_threshold = d_absolute_residual_tolerance_threshold.getValue();
        const auto & newton_iterations = d_newton_iterations.getValue();
        const auto & print_log = f_printLog.getValue();
        auto info = MessageDispatcher::info(Message::Runtime, ComponentInfo::SPtr(new ComponentInfo(this->getClassName())), SOFA_FILE_INFO);

        // We are (normally) at the beginning of the time step, let's specify that we haven't analyze the pattern of the system matrix yet
        if (pattern_strategy == PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP) {
            p_has_already_analyzed_the_pattern = false; // This will allow the matrix to be analyzed after the next assembly
        }

        // Right hand side term (internal + external forces)
        auto f_id = sofa::core::MultiVecDerivId(sofa::core::VecDerivId::force());
        vop.v_clear(f_id);

        // Incremental displacement of one iteration (not allocated by default by the mechanical objects, unlike x, v, f and df)
        auto dx_id = sofa::core::MultiVecDerivId(sofa::core::VecDerivId::dx());
        vop.v_realloc(dx_id, false /* interactionForceField */, false /* propagate [to mapped MO] */);
        vop.v_clear(dx_id);

        // Total displacement increment since the beginning
        vop.v_realloc(p_U_id, false /* interactionForceField */, false /* propagate [to mapped MO] */);
        vop.v_clear(p_U_id);

        // Set implicit param to true to trigger nonlinear stiffness matrix recomputation
        mop->setImplicit(true);

        if (print_log) {
            info << "======= Starting static ODE solver =======\n";
            info << "Time step                : " << this->getTime() << "\n";
            info << "Context                  : " << dynamic_cast<const sofa::simulation::Node *>(context)->getPathName() << "\n";
            info << "Max iterations           : " << newton_iterations << "\n";
            info << "Residual tolerance (abs) : " << absolute_residual_tolerance_threshold << "\n";
            info << "Residual tolerance (rel) : " << residual_tolerance_threshold << "\n";
            info << "Correction tolerance     : " << correction_tolerance_threshold << "\n";
            info << "Linear solver            : " << l_linear_solver->getPathName() << "\n\n";
        }

        // Local variables used for the iterations
        unsigned n_it=0;
        double dx_squared_norm, du_squared_norm, R_squared_norm = 0;
        const auto squared_residual_threshold = residual_tolerance_threshold*residual_tolerance_threshold;
        const auto squared_correction_threshold = correction_tolerance_threshold*correction_tolerance_threshold;
        const auto squared_absolute_residual_tolerance_threshold = absolute_residual_tolerance_threshold*absolute_residual_tolerance_threshold;
        bool converged = false, diverged = false;
        steady_clock::time_point t;

        // Resize vectors containing the newton residual norms
        p_squared_residuals.clear();
        p_squared_residuals.reserve(newton_iterations);

        // Resize vectors containing the times took to compute the newton iterations
        p_times.clear();
        p_times.reserve(newton_iterations);

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

        p_DX.reset(linear_solver->create_new_vector(n));
        p_DX->clear();

        p_F.reset(linear_solver->create_new_vector(n));
        p_F->clear();


        // ###########################################################################
        // #                             First residual                              #
        // ###########################################################################
        // # Before starting any newton iterations, we first need to compute         #
        // # the residual with the updated right-hand side (the new load increment)  #
        // ###########################################################################

        // Step 1   Assemble the force vector
        sofa::helper::AdvancedTimer::stepBegin("ComputeForce");
        this->assemble_rhs_vector(mechanical_parameters, accessor, f_id, p_F.get());
        sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

        // Step 2   Compute the initial residual
        R_squared_norm = SofaCaribou::Algebra::dot(p_F.get(), p_F.get());
        p_squared_initial_residual = R_squared_norm;

        if (absolute_residual_tolerance_threshold > 0 && R_squared_norm <= squared_absolute_residual_tolerance_threshold) {
            converged = true;
            if (print_log) {
                info << "The ODE has already reached an equilibrium state."
                     << std::scientific
                     << "The residual's ratio |R| is " << std::setw(12) << sqrt(R_squared_norm)
                     << " (criterion is " << std::setw(12) << absolute_residual_tolerance_threshold << ") \n"
                     << std::defaultfloat;
            }
        }

        // ###########################################################################
        // #                          Newton iterations                              #
        // ###########################################################################

        while (not converged and n_it < newton_iterations) {
            sofa::helper::ScopedAdvancedTimer step_timer ("NewtonStep");
            t = steady_clock::now();

            // Part 1. Assemble the system matrix.
            {
                sofa::helper::ScopedAdvancedTimer _t_("MBKBuild");
                p_A->clear();
                this->assemble_system_matrix(mechanical_parameters, accessor, p_A.get());
                linear_solver->set_system_matrix(p_A.get());
            }

            // Part 2. Analyze the pattern of the matrix in order to compute a permutation matrix.
            {
                // Let's see if we should (re)-analyze the pattern of the system matrix
                if (
                        pattern_strategy != PatternAnalysisStrategy::NEVER and (
                                pattern_strategy == PatternAnalysisStrategy::ALWAYS or
                                (
                                        (pattern_strategy == PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP or pattern_strategy == PatternAnalysisStrategy::BEGINNING_OF_THE_SIMULATION)
                                        and not p_has_already_analyzed_the_pattern
                                )
                        )
                        ) {

                    sofa::helper::ScopedAdvancedTimer _t_("MBKAnalyze");

                    if (not linear_solver->analyze_pattern()) {
                        info << "[DIVERGED] Failed to analyze the pattern of the system matrix.";
                        diverged = true;
                        break;
                    }

                    p_has_already_analyzed_the_pattern = true;
                }
            }

            // Part 3. Factorize the matrix.
            {
                sofa::helper::ScopedAdvancedTimer _t_("MBKFactorize");
                if (not linear_solver->factorize()) {
                    info << "[DIVERGED] Failed to factorize the system matrix.";
                    diverged = true;
                    break;
                }
            }

            // Part 4. Solve the unknown increment.
            {
                sofa::helper::ScopedAdvancedTimer _t_("MBKSolve");
                if (not linear_solver->solve(p_F.get(), p_DX.get())) {
                    info << "[DIVERGED] The linear solver failed to solve the unknown increment.";
                    diverged = true;
                    break;
                }
            }

            // Part 5. Propagating the solution increment and update geometry.
            {
                sofa::helper::ScopedAdvancedTimer _t_("PropagateDx");
                this->propagate_solution_increment(mechanical_parameters, accessor, p_DX.get(), x_id, v_id, dx_id);
            }

            // The next two parts are only necessary when doing more than one Newton iteration
            if (newton_iterations > 1) {
                // Part 6. Update the force vector.
                sofa::helper::AdvancedTimer::stepBegin("UpdateForce");
                p_F->clear();
                this->assemble_rhs_vector(mechanical_parameters, accessor, f_id, p_F.get());
                sofa::helper::AdvancedTimer::stepEnd("UpdateForce");

                // Part 7. Compute the updated force residual.
                sofa::helper::AdvancedTimer::stepBegin("UpdateResidual");
                R_squared_norm = SofaCaribou::Algebra::dot(p_F.get(), p_F.get());
                sofa::helper::AdvancedTimer::stepEnd("UpdateResidual");
            }

            // Part 8. Compute the updated displacement residual.
            sofa::helper::AdvancedTimer::stepBegin("UpdateU");
            vop.v_peq(p_U_id, dx_id); // U += dx
            vop.v_dot(dx_id, dx_id);  // dx.dot(dx)
            dx_squared_norm = vop.finish();

            vop.v_dot(p_U_id, p_U_id); // U.dot(U)
            du_squared_norm = vop.finish();
            sofa::helper::AdvancedTimer::stepEnd("UpdateU");

            // Part 9. Stop timers and print step information.
            auto iteration_time = duration_cast<nanoseconds>(steady_clock::now() - t).count();
            p_times.emplace_back(static_cast<UNSIGNED_INTEGER_TYPE>(iteration_time));

            p_squared_residuals.emplace_back(R_squared_norm);

            // We completed one iteration, increment the counter
            n_it++;

            if( print_log ) {
                info << "Newton iteration #" << std::left << std::setw(5)  << n_it
                     << std::scientific
                     << "  |R| = "   << std::setw(12) << sqrt(R_squared_norm)
                     << "  |R|/|R0| = "   << std::setw(12) << sqrt(R_squared_norm  / p_squared_residuals[0])
                     << "  |du| / |U| = " << std::setw(12) << sqrt(dx_squared_norm / du_squared_norm)
                     << std::defaultfloat;
                info << "  Time = " << iteration_time/1000/1000 << " ms";
                if (linear_solver->is_iterative()) {
                    info << "  # of linear solver iterations = " << linear_solver->squared_residuals().size();
                }
                info << "\n";
            }

            if (std::isnan(R_squared_norm) or std::isnan(dx_squared_norm) or du_squared_norm < EPSILON) {
                diverged = true;
                if (print_log) {
                    info << "[DIVERGED]";
                    if (std::isnan(R_squared_norm)) {
                        info << " The residual's ratio |R| is NaN.";
                    }
                    if (std::isnan(dx_squared_norm)) {
                        info << " The correction's ratio |du| is NaN.";
                    }
                    if (du_squared_norm < EPSILON) {
                        info << " The correction's ratio |du|/|U| is NaN (|U| is zero).";
                    }
                    info << "\n";
                }
                break;
            }

            // Part 10. Check for convergence.
            if (correction_tolerance_threshold > 0 and dx_squared_norm < squared_correction_threshold*du_squared_norm) {
                converged = true;
                if (print_log) {
                    info  << "[CONVERGED] The correction's ratio |du|/|U| = " << sqrt(dx_squared_norm/du_squared_norm) << " is smaller than the threshold of " << correction_tolerance_threshold << ".\n";
                }
                break;
            }

            if (residual_tolerance_threshold > 0 and R_squared_norm < squared_residual_threshold*p_squared_residuals[0]) {
                converged = true;
                if (print_log) {
                    info << "[CONVERGED] The residual's ratio |R|/|R0| = " << sqrt(R_squared_norm/p_squared_residuals[0]) << " is smaller than the threshold of " << residual_tolerance_threshold << ".\n";
                }
                break;
            }

            if (absolute_residual_tolerance_threshold > 0 and R_squared_norm < squared_absolute_residual_tolerance_threshold) {
                converged = true;
                if (print_log) {
                    info << "[CONVERGED] The residual's ratio |R| = " << sqrt(R_squared_norm) << " is smaller than the threshold of " << absolute_residual_tolerance_threshold << ".\n";
                }
                break;
            }

            // Clear up the solution
            vop.v_clear(dx_id);
        } // End while (not converged and not diverged and n_it < newton_iterations)

        n_it--; // Reset to the actual index of the last iteration completed

        if (not converged and not diverged and n_it == (newton_iterations-1)) {
            if (print_log) {
                info << "[DIVERGED] The number of Newton iterations reached the maximum of " << newton_iterations << " iterations" << ".\n";
            }
        }

        d_converged.setValue(converged);

        sofa::helper::AdvancedTimer::valSet("has_converged", converged ? 1 : 0);
        sofa::helper::AdvancedTimer::valSet("nb_iterations", n_it+1);
    }


/**
 Copy pasted from https://github.com/jnbrunet/caribou/blob/master/src/SofaCaribou/Ode/StaticODESolver.cpp
 */
// Assemble F in A [dx] = F
    void NewtonRaphson::assemble_rhs_vector(const sofa::core::MechanicalParams &mechanical_parameters,
                                              const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                              sofa::core::MultiVecDerivId & f_id,
                                              sofa::defaulttype::BaseVector *f) {
        // 1. Clear the force vector (F := 0)
        MechanicalResetForceVisitor(&mechanical_parameters, f_id,
                                    false /* onlyMapped */)
                .execute(this->getContext());

        // 2. Go down in the current context tree calling `addForce` on every force field components,
        //    then go up from the leaves calling `applyJT` on every mechanical mappings
        MechanicalComputeForceVisitor(&mechanical_parameters, f_id,
                                      true /* accumulate (to mapped node) */, true /*neglectingCompliance*/)
                .execute(this->getContext());

        // 3. Calls the "projectResponse" method of every `BaseProjectiveConstraintSet` objects found in the current
        //    context tree. For example, the `FixedConstraints` component will set entries of fixed nodes to zero.
        MechanicalApplyConstraintsVisitor(&mechanical_parameters, f_id,
                                          nullptr /*W (also project the given compliance matrix) */)
                .execute(this->getContext());

        // 4. Copy force vectors from every top level (unmapped) mechanical objects into the given system vector f
        MechanicalMultiVectorToBaseVectorVisitor(&mechanical_parameters, f_id /* source */, f /* destination */, &matrix_accessor)
                .execute(this->getContext());
    }

// Assemble A in A [dx] = F
    void NewtonRaphson::assemble_system_matrix(const sofa::core::MechanicalParams &mechanical_parameters,
                                                 sofa::component::linearsolver::DefaultMultiMatrixAccessor & matrix_accessor,
                                                 sofa::defaulttype::BaseMatrix *A) {
        // Step 1. Building stage:
        //         Here we go down on the current context sub-graph and call :
        //           1. ff->addMBKToMatrix(&K) for every force field "ff" found.
        //           2. pc->applyConstraint(&K) for every BaseProjectiveConstraintSet "pc" found.
        //         If a mechanical mapping "m" is found during the traversal, and m->areMatricesMapped() is false, the
        //         traversal stops in the subgraph of the mapping.
        matrix_accessor.setGlobalMatrix(A);
        sofa::core::MechanicalParams m_params (mechanical_parameters);
        m_params.setKFactor(-1.0);
        Timer::stepBegin("AssembleGlobalMatrix");
        SofaCaribou::visitor::AssembleGlobalMatrix(&m_params, &matrix_accessor).execute(this->getContext());
        Timer::stepEnd("AssembleGlobalMatrix");

        Timer::stepBegin("ConstrainGlobalMatrix");
        SofaCaribou::visitor::ConstrainGlobalMatrix(&m_params, &matrix_accessor).execute(this->getContext());
        Timer::stepEnd("ConstrainGlobalMatrix");

        // Step 2. Mechanical mappings
        //         In case we have mapped matrices, which is, system matrix of a slave mechanical object, accumulate its
        //         contribution to the global system matrix with:
        //           [A]ij += Jt * [A']ij * J
        //         where A is the master mechanical object's matrix, A' is the slave mechanical object matrix and J=m.getJ()
        //         is the mapping relation between the slave and its master.
        Timer::stepBegin("MappedMatrices");
        matrix_accessor.computeGlobalMatrix();
        Timer::stepEnd("MappedMatrices");

        // Step 3. Convert the system matrix to a compressed sparse matrix
        Timer::stepBegin("ConvertToSparse");
        A->compress();
        Timer::stepEnd("ConvertToSparse");

    }

// Propagate Dx that was previously solved in A [dx] = F
    void NewtonRaphson::propagate_solution_increment(const sofa::core::MechanicalParams &mechanical_parameters,
                                                       const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                                       const sofa::defaulttype::BaseVector * dx,
                                                       sofa::core::MultiVecCoordId & x_id,
                                                       sofa::core::MultiVecDerivId & /* v_id */,
                                                       sofa::core::MultiVecDerivId & dx_id) {

        // 1. Copy vectors from the global system vector into every top level (unmapped) mechanical objects.
        MechanicalMultiVectorFromBaseVectorVisitor(&mechanical_parameters, dx_id, dx, &matrix_accessor).execute(this->getContext());

        // 2. x += dx
        MechanicalVOpVisitor(&mechanical_parameters, x_id, x_id, dx_id).execute(this->getContext());

        // 3. Calls "solveConstraint" method of every ConstraintSolver objects found in the current context tree.
        sofa::core::ConstraintParams constraint_parameters = mechanical_parameters;
        constraint_parameters.setOrder(sofa::core::ConstraintParams::POS);

        using Direction = sofa::core::objectmodel::BaseContext::SearchDirection;
        auto constraint_solvers = this->getContext()->getObjects<sofa::core::behavior::ConstraintSolver>(Direction::Local);
        for (auto * solver : constraint_solvers) {
            solver->solveConstraint(&constraint_parameters, x_id);
        }

        // 4. Propagate positions to mapped mechanical objects, for example, identity mappings, barycentric mappings, etc.
        //    This will call the methods apply and applyJ on every mechanical mappings.
        MechanicalPropagateOnlyPositionAndVelocityVisitor(&mechanical_parameters).execute(this->getContext());
    }

    int NewtonRaphsonClass = sofa::core::RegisterObject("Newton Raphson").add< NewtonRaphson >();
} // namespace DeepPhysicsSofa::ode