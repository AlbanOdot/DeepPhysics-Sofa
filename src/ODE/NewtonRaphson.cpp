
/**
 * This class is hugely inspired by https://github.com/jnbrunet/caribou/blob/master/src/SofaCaribou/Ode/NewtonRaphsonSolver.cpp
 */

#include <src/ODE/NewtonRaphson.h>

#include <iomanip>
#include <chrono>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/VectorOperations.h>
        DISABLE_ALL_WARNINGS_BEGIN

#include <SofaCaribou/Solver/LinearSolver.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Size = int; }
#endif

namespace DeepPhysicsSofa::ode {

    using sofa::core::ExecParams;
    using sofa::core::MultiVecCoordId;
    using sofa::core::MultiVecDerivId;

    NewtonRaphson::NewtonRaphson()
            : d_newton_iterations(initData(&d_newton_iterations,
                                           (unsigned) 1,
                                           "newton_iterations",
                                           "Number of newton iterations between each load increments (normally, one load increment per simulation time-step."))
            , d_correction_tolerance_threshold(initData(&d_correction_tolerance_threshold,
                                                        (double) 1e-5,
                                                        "correction_tolerance_threshold",
                                                        "Relative convergence criterion: The newton iterations will stop when the norm of correction |du| reach this threshold."))
            , d_residual_tolerance_threshold( initData(&d_residual_tolerance_threshold,
                                                       (double) 1e-5,
                                                       "residual_tolerance_threshold",
                                                       "Relative convergence criterion: The newton iterations will stop when the ratio between norm of the residual "
                                                       "R_k = |f_k - K(u_k)| at iteration k over R_0 is lower than this threshold. Use a negative value to "
                                                       "disable this criterion."))
            , d_absolute_residual_tolerance_threshold( initData(&d_absolute_residual_tolerance_threshold,
                                                                (double) 1e-15,
                                                                "absolute_residual_tolerance_threshold",
                                                                "Absolute convergence criterion: The newton iterations will stop when the norm of the residual "
                                                                "R_k = |f_k - K(u_k)| at iteration k is lower than this threshold. Use a negative value to "
                                                                "disable this criterion."))
            , d_pattern_analysis_strategy(initData(&d_pattern_analysis_strategy,
                                                   "pattern_analysis_strategy",
                                                   "Define when the pattern of the system matrix should be analyzed to extract a permutation matrix. If the sparsity and"
                                                   "location of the coefficients of the system matrix doesn't change much during the simulation, then this analysis can"
                                                   "be avoided altogether, or computed only one time at the beginning of the simulation. Else, it can be done at the "
                                                   "beginning of the time step, or even at each reformation of the system matrix if necessary. The default is to "
                                                   "analyze the pattern at each time step."))
            , l_linear_solver(initLink(
                    "linear_solver",
                    "Linear solver used for the resolution of the system."))
            , d_converged(initData(&d_converged,
                                   false,
                                   "converged",
                                   "Whether or not the last call to solve converged",
                                   true /*is_displayed_in_gui*/,
                                   true /*is_read_only*/))
    {
        d_pattern_analysis_strategy.setValue(sofa::helper::OptionsGroup(std::vector < std::string > {
                "NEVER", "BEGINNING_OF_THE_SIMULATION", "BEGINNING_OF_THE_TIME_STEP", "ALWAYS"
        }));

        // Select the default value
        set_pattern_analysis_strategy(PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP);
    }
    void NewtonRaphson::init() {
        p_has_already_analyzed_the_pattern = false;

        if (not has_valid_linear_solver()) {
            // No linear solver specified, let's try to find one in the current node
            auto solvers = this->getContext()->template getObjects<sofa::core::behavior::LinearSolver>(sofa::core::objectmodel::BaseContext::Local);
            std::vector<sofa::core::behavior::LinearSolver *> sofa_linear_solvers;
            std::vector<sofa::core::behavior::LinearSolver *> caribou_linear_solvers;

            // If we have many linear solvers in the current node, let's try to classify them in SOFA vs Caribou
            for(auto * solver : solvers) {
                auto caribou_solver = dynamic_cast<SofaCaribou::solver::LinearSolver *> (solver);
                if (caribou_solver) {
                    caribou_linear_solvers.push_back(solver);
                } else {
                    sofa_linear_solvers.push_back(solver);
                }
            }


            if (caribou_linear_solvers.empty()) {
                if (sofa_linear_solvers.empty()) {
                    // No caribou and no SOFA linear solvers were found.
                    msg_error() << "No compatible linear solvers were found in the current context. The '"
                                << l_linear_solver.getName()
                                << "' attribute can be use to specify the path towards a linear solver.";
                } else {
                    // No caribou linear solver found, but SOFA's linear solver were found. Maybe a user mistake?
                    msg_error() << sofa_linear_solvers.size()
                                << " linear solver were found, none of which are compatible with this ODE solver. The '"
                                << l_linear_solver.getName()
                                << "' attribute can be use to specify the path towards a compatible linear solver.";
                }
            } else if (caribou_linear_solvers.size() == 1) {
                // If we have only one Caribou linear solver, let's take it
                l_linear_solver.set(caribou_linear_solvers[0]);
                msg_info() << "Automatically found the linear solver '"
                           << l_linear_solver->getPathName()
                           << "' from the current context. If another one was expected, use the '"
                           << l_linear_solver.getName()
                           << "' attribute.";
            } else {
                // Multiple Caribou linear solver were found... Let's take the first one and notify the user that it might
                // not be the good one.
                l_linear_solver.set(caribou_linear_solvers[0]);
                msg_warning() << "Multiple compatible linear solvers were found in the current context. The first one ("
                              << l_linear_solver->getPathName() << ") will be used. If another one was expected, or to "<<
                              "remove of this warning, use the '" << l_linear_solver.getName() << "' attribute.";
            }

        }
    }

    void NewtonRaphson::reset() {
        p_has_already_analyzed_the_pattern = false;
    }

    bool NewtonRaphson::has_valid_linear_solver() const {
        return (
                l_linear_solver.get() != nullptr and
                dynamic_cast<SofaCaribou::solver::LinearSolver *> (l_linear_solver.get()) != nullptr
        );
    }

    auto NewtonRaphson::pattern_analysis_strategy() const -> NewtonRaphson::PatternAnalysisStrategy {
        const auto v = static_cast<PatternAnalysisStrategy>(d_pattern_analysis_strategy.getValue().getSelectedId());
        switch (v) {
            case PatternAnalysisStrategy::ALWAYS:
            case PatternAnalysisStrategy::BEGINNING_OF_THE_SIMULATION:
            case PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP:
            case PatternAnalysisStrategy::NEVER:
                return v;
        }

        // Default value
        return NewtonRaphson::PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP;
    }

    void NewtonRaphson::set_pattern_analysis_strategy(const NewtonRaphson::PatternAnalysisStrategy & strategy) {
        using namespace sofa::helper;
        auto pattern_analysis_strategy = WriteOnlyAccessor<Data<OptionsGroup>>(d_pattern_analysis_strategy);
        pattern_analysis_strategy->setSelectedItem(static_cast<unsigned int> (strategy));
    }

} // namespace DeepPhysicsSofa::ode