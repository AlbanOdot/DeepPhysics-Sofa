//
// Created by alban on 14/03/2022.
//

#ifndef DEEPPHYSICSSOFA_KSOLVER_H
#define DEEPPHYSICSSOFA_KSOLVER_H
#include <SofaCaribou/Ode/StaticODESolver.h>
namespace DeepPhysicsSofa {
    namespace ode {

        class KSolver : public SofaCaribou::ode::StaticODESolver {
        public:
            SOFA_CLASS(KSolver, SofaCaribou::ode::StaticODESolver);
            KSolver();
            template<typename T>
            using Data = sofa::core::objectmodel::Data<T>;

        protected:
            void solve(const sofa::core::ExecParams *params, SReal dt, sofa::core::MultiVecCoordId x_id,
                       sofa::core::MultiVecDerivId v_id) override;

            // Value of the residual from the prediction state
            Data<double> d_prediction_residual;
        };

    }
}

#endif //DEEPPHYSICSSOFA_KSOLVER_H
