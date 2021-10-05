#pragma once

#include <SofaBoundaryCondition/ConstantForceField.h>
namespace DeepPhysicsSofa::forcefield {


    template<class DataTypes>
    class FreeConstantForceField : public sofa::component::forcefield::ConstantForceField<DataTypes> {
    public:
        SOFA_CLASS(SOFA_TEMPLATE(FreeConstantForceField, DataTypes), SOFA_TEMPLATE(
                sofa::component::forcefield::ConstantForceField, DataTypes)

        );
    public:
        typedef sofa::core::behavior::ForceField <DataTypes> Inherit;
        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;
        typedef sofa::type::vector<unsigned int> VecIndex;
        typedef sofa::core::objectmodel::Data <VecCoord> DataVecCoord;
        typedef sofa::core::objectmodel::Data <VecDeriv> DataVecDeriv;

        typedef sofa::component::topology::PointSubsetData <VecIndex> SetIndex;

        /// Set a force to a given particle
        void setForce(unsigned i, const Deriv &f);

    protected:
        FreeConstantForceField();

    };

}
