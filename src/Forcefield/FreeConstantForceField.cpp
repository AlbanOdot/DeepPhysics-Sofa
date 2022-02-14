#include <src/Forcefield/FreeConstantForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <numeric>

namespace DeepPhysicsSofa::forcefield {
    using namespace sofa::core;
    using objectmodel::BaseObjectDescription ;

    template<class DataTypes>
    FreeConstantForceField<DataTypes>::FreeConstantForceField():sofa::component::forcefield::ConstantForceField<DataTypes> ()
    {
        this->d_showArrowSize.setGroup("Visualization");
        this->d_color.setGroup("Visualization");
    }

    template<class DataTypes>
    void FreeConstantForceField<DataTypes>::init()
    {

        if (this->l_topology.empty())
        {
            msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
            this->l_topology.set(this->getContext()->getMeshTopologyLink());
        }

        // temprory pointer to topology
        topology::BaseMeshTopology* _topology = this->l_topology.get();

        if (_topology)
        {
            msg_info() << "Topology path used: '" << this->l_topology.getLinkedPath() << "'";

            // Initialize functions and parameters for topology data and handler
            this->d_indices.createTopologyHandler(_topology);

            m_systemSize = _topology->getNbPoints();
        }
        else
        {
            msg_info() << "No topology component found at path: " << this->l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name;
            behavior::BaseMechanicalState* state = this->getContext()->getMechanicalState();
            m_systemSize = state->getSize();
        }


        const VecIndex & indices = this->d_indices.getValue();
        size_t indicesSize = indices.size();

        if (this->d_indices.isSet() && indicesSize!=0)
        {
            // check size of vector indices
            if( indicesSize > m_systemSize )
            {
                msg_error() << "Size mismatch: indices > system size";

                return;
            }
            // check each indice of the vector
            for(size_t i=0; i<indicesSize; i++)
            {
                if( indices[i] > m_systemSize )
                {
                    msg_error() << "Indices incorrect: indice["<< i <<"] = "<< indices[i] <<" exceeds system size";

                    return;
                }
            }
        }
        else
        {
            // initialize with all indices
            VecIndex& indicesEdit = *this->d_indices.beginEdit();
            indicesEdit.clear();
            indicesEdit.resize(m_systemSize);
            std::iota (std::begin(indicesEdit), std::end(indicesEdit), 0);
            this->d_indices.endEdit();
        }

        if (this->d_forces.isSet())
        {
            const VecDeriv &forces = this->d_forces.getValue();
            if( !checkForces(forces) )
            {
                msg_error() << " Invalid given vector forces";

                return;
            }
            msg_info() << "Input vector forces is used for initialization";
        }
        else if (this->d_force.isSet())
        {
            const Deriv &force = this->d_force.getValue();
            if( checkForce(force) )
            {
                computeForceFromSingleForce();
            }
            else
            {
                msg_error() << " Invalid given force";

                return;
            }
            msg_info() << "Input force is used for initialization";
        }

        // init from ForceField
        Inherit::init();

        // add to tracker
        this->trackInternalData(this->d_indices);
        this->trackInternalData(this->d_forces);
        this->trackInternalData(this->d_force);

        // if all init passes, component is valid

    }

    template<class DataTypes>
    void FreeConstantForceField<DataTypes>::doUpdateInternal()
    {
        if (this->hasDataChanged(this->d_indices))
        {
            msg_info() << "doUpdateInternal: data indices has changed";

            const VecIndex & indices = this->d_indices.getValue();
            size_t indicesSize = indices.size();

            // check size of vector indices
            if( indicesSize > m_systemSize )
            {
                msg_error() << "Size mismatch: indices > system size";

                return;
            }
            else if( indicesSize==0 )
            {
                msg_warning() << "Size of vector indices is zero";
            }

            // check each indice of the vector
            for(size_t i=0; i<indicesSize; i++)
            {
                if( indices[i] > m_systemSize )
                {
                    msg_error() << "Indices incorrect: indice["<< i <<"] = "<< indices[i] <<" exceeds system size";

                    return;
                }
            }
        }

        if (this->hasDataChanged(this->d_forces))
        {
            msg_info() << "doUpdateInternal: data forces has changed";

            const VecDeriv &forces = this->d_forces.getValue();
            if( checkForces(forces) )
            {
                computeForceFromForceVector();

            }
            else
            {
                msg_error() << " Invalid given vector forces";

                return;
            }
        }

        if (this->hasDataChanged(this->d_force))
        {
            msg_info() << "doUpdateInternal: data force has changed";

            const Deriv &force = this->d_force.getValue();
            if( checkForce(force) )
            {
                computeForceFromSingleForce();

            }
            else
            {
                msg_error() << " Invalid given force";

                return;
            }
        }
    }



    template<class DataTypes>
    void FreeConstantForceField<DataTypes>::computeForceFromSingleForce()
    {
        const Deriv& singleForce = this->d_force.getValue();
        const VecIndex & indices = this->d_indices.getValue();
        VecDeriv& forces = *this->d_forces.beginEdit();
        size_t indicesSize = indices.size();
        forces.clear();
        forces.resize(indicesSize);

        for(size_t i=0; i<indicesSize; i++)
        {
            forces[i] = singleForce;
        }

        this->d_forces.endEdit();
    }

    template<class DataTypes>
    void FreeConstantForceField<DataTypes>::computeForceFromForceVector()
    {
        const VecDeriv& forces = this->d_forces.getValue();
        const size_t indicesSize = this->d_indices.getValue().size();
        Deriv& totalForce = *this->d_totalForce.beginEdit();
        totalForce.clear();
        if( indicesSize!=forces.size() )
        {
            msg_error() << "Impossible to use the vector forces since its size mismatches with indices size";
            return;
        }
        else
        {
            for(size_t i=0; i<indicesSize; i++)
            {
                totalForce += forces[i];
            }
        }

        this->d_totalForce.endEdit();
    }

    template <class DataTypes>
    void FreeConstantForceField<DataTypes>::setForce(unsigned i, const Deriv& force)
    {
        VecIndex& indices = *this->d_indices.beginEdit();
        VecDeriv& f = *this->d_forces.beginEdit();
        indices.push_back(i);
        f.push_back( force );
        this->d_indices.endEdit();
        this->d_forces.endEdit();
    }

    template<class DataTypes>
    bool FreeConstantForceField<DataTypes>::checkForce(const Deriv& force)
    {
        size_t size = Deriv::spatial_dimensions;

        for (size_t i=0; i<size; i++)
        {
            if( std::isnan(force[i]) )
            {
                return false;
            }
        }
        return true;
    }

    template<class DataTypes>
    bool FreeConstantForceField<DataTypes>::checkForces(const VecDeriv& forces)
    {
        for (auto&& i : forces)
        {
            if(! checkForce(i) )
            {
                return false;
            }
        }
        return true;
    }

    using namespace sofa::defaulttype;

    int FreeConstantForceFieldClass = sofa::core::RegisterObject("Constant forces applied to given degrees of freedom")
                                              .add< FreeConstantForceField<Vec3Types> >()
                                              .add< FreeConstantForceField<Vec2Types> >()
                                              .add< FreeConstantForceField<Vec1Types> >()
                                              .add< FreeConstantForceField<Vec6Types> >()
                                              .add< FreeConstantForceField<Rigid3Types> >()
                                              .add< FreeConstantForceField<Rigid2Types> >()

    ;

}