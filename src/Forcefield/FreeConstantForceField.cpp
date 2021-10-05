#include <src/Forcefield/FreeConstantForceField.h>
#include <sofa/core/ObjectFactory.h>

namespace DeepPhysicsSofa::forcefield {
    using namespace sofa::core;
    using objectmodel::BaseObjectDescription;

    template<class DataTypes>
    FreeConstantForceField<DataTypes>::FreeConstantForceField()
            :sofa::component::forcefield::ConstantForceField<DataTypes>() {}

    template<class DataTypes>
    void FreeConstantForceField<DataTypes>::setForce(unsigned i, const Deriv &force) {
        VecIndex &indices = *this->d_indices.beginEdit();
        VecDeriv &f = *this->d_forces.beginEdit();
        indices.push_back(i);
        f.push_back(force);
        this->d_indices.endEdit();
        this->d_forces.endEdit();
    }

    using namespace sofa::defaulttype;

    int FreeConstantForceFieldClass = sofa::core::RegisterObject("Constant forces applied to given degrees of freedom")
                                              .add < FreeConstantForceField < Vec3Types > > ()
                                              .add < FreeConstantForceField < Vec2Types > > ()
                                              .add < FreeConstantForceField < Vec1Types > > ()
                                              .add < FreeConstantForceField < Vec6Types > > ()
                                              .add < FreeConstantForceField < Rigid3Types > > ()
                                              .add < FreeConstantForceField < Rigid2Types > > ();
}