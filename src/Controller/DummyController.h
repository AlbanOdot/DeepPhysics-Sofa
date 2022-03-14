//
// Created by alban on 11/03/2022.
//

#ifndef DEEPPHYSICSSOFA_DUMMYCONTROLLER_H
#define DEEPPHYSICSSOFA_DUMMYCONTROLLER_H
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/simulation/AnimateBeginEvent.h>
class DummyController : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(DummyController, BaseObject);
    using DataTypes = sofa::defaulttype::Vec3Types;

    sofa::core::objectmodel::SingleLink< DummyController,
                                         sofa::core::behavior::MechanicalState<DataTypes>,
                                         sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
                                         sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK> l_mstate;

    DummyController()
            : l_mstate(initLink("MO", "Mechanical state to reorder"))
    {
        this->f_listening.setValue(true);
    }

    void handleEvent(sofa::core::objectmodel::Event* event)
    {
        const auto& moValue = sofa::helper::getWriteAccessor(*l_mstate.get()->write(sofa::core::VecCoordId::position()));
    }
};

#endif //DEEPPHYSICSSOFA_DUMMYCONTROLLER_H
