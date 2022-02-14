#pragma once

#include <SofaBoundaryCondition/ConstantForceField.h>
#include <sofa/core/topology/TopologySubsetIndices.h>
namespace DeepPhysicsSofa::forcefield {


    template<class DataTypes>
    class FreeConstantForceField : public sofa::component::forcefield::ConstantForceField<DataTypes> {
    public:
        SOFA_CLASS(SOFA_TEMPLATE(FreeConstantForceField, DataTypes), SOFA_TEMPLATE(sofa::component::forcefield::ConstantForceField, DataTypes));
        typedef sofa::core::behavior::ForceField<DataTypes> Inherit;
        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;
        typedef std::vector<unsigned int> VecIndex;
        typedef sofa::core::objectmodel::Data<VecCoord> DataVecCoord;
        typedef sofa::core::objectmodel::Data<VecDeriv> DataVecDeriv;

        typedef sofa::core::topology::TopologySubsetIndices SetIndex;
    public:

        /// Init function
        void init() override;

        /// Update data and internal vectors
        void doUpdateInternal() override;

        /// Set a force to a given particle
        void setForce( unsigned i, const Deriv& f );


        using Inherit::addAlias ;
        using Inherit::addKToMatrix;


    protected:
        FreeConstantForceField();

        /// Functions computing and updating the constant force vector
        void computeForceFromSingleForce();
        void computeForceFromForceVector();

        /// Functions checking inputs before update
        bool checkForce(const Deriv&  force);
        bool checkForces(const VecDeriv& forces);

        /// Save system size for update of indices (doUpdateInternal)
        size_t m_systemSize;
    };


}
