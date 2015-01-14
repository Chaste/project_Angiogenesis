/*
 * CaCaVascularNetworkNode.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#include "CaVascularNetworkNode.hpp"
#include "CaVessel.hpp"
#include <cassert>


template<unsigned SPATIAL_DIM>
CaVascularNetworkNode<SPATIAL_DIM>::CaVascularNetworkNode():mLocation(),mPressure(0.0),pAdjoiningVessels(),mIsInputNode(false),mIsOutputNode(false)
{

}

template<unsigned SPATIAL_DIM>
CaVascularNetworkNode<SPATIAL_DIM>::~CaVascularNetworkNode()
{

}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > CaVascularNetworkNode<SPATIAL_DIM>::shared()
{
    return this->shared_from_this();
}

template<unsigned SPATIAL_DIM>
ChastePoint<SPATIAL_DIM> CaVascularNetworkNode<SPATIAL_DIM>::GetLocation()
{
    return mLocation;
}

template<unsigned SPATIAL_DIM>
double CaVascularNetworkNode<SPATIAL_DIM>::GetPressure() const
{
    return mPressure;
}

template<unsigned SPATIAL_DIM>
int CaVascularNetworkNode<SPATIAL_DIM>::GetNumberOfAdjoiningVessels() const
{
    return pAdjoiningVessels.size();
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > CaVascularNetworkNode<SPATIAL_DIM>::GetAdjoiningVessel(int i)
{
    return pAdjoiningVessels[i].lock(); // lock() converts weak pointer (stored here) to shared pointer
}

template<unsigned SPATIAL_DIM>
bool CaVascularNetworkNode<SPATIAL_DIM>::IsInputNode()
{
    return mIsInputNode;
}

template<unsigned SPATIAL_DIM>
bool CaVascularNetworkNode<SPATIAL_DIM>::IsOutputNode()
{
    return mIsOutputNode;
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::SetLocation(ChastePoint<SPATIAL_DIM> loc)
{
    mLocation = loc;
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::SetPressure(double pressure)
{
    mPressure = pressure;
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::AddAdjoiningVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{

    // same vessel may be adjoint to a single node twice but that node must be both node1 and node2 of the doubly adjoint vessel

    int numberOfTimeVesselIsAlreadyAttachedToNode = 0;

    for(int i = 0; i < pAdjoiningVessels.size(); i++)
    {
        if (pAdjoiningVessels[i].lock() == vessel)
        {
            numberOfTimeVesselIsAlreadyAttachedToNode++;
        }
    }

    assert(numberOfTimeVesselIsAlreadyAttachedToNode <= 2);

    if (numberOfTimeVesselIsAlreadyAttachedToNode == 0)
    {
        pAdjoiningVessels.push_back(boost::weak_ptr<CaVessel<SPATIAL_DIM> >(vessel));
    }
    else if (numberOfTimeVesselIsAlreadyAttachedToNode == 1 && vessel->GetNode1() == shared() && vessel->GetNode2() == shared())
    {
        pAdjoiningVessels.push_back(boost::weak_ptr<CaVessel<SPATIAL_DIM> >(vessel));
    }
    else
    {
        // vessel is already attached to node correctly (once or twice)
    }

    /*
        todo should check that vessel being adjoined to node has not got an actively migrating tip located at this node
        todo should perhaps check that the vessel exists at the location where the node is located
        - i.e. that the vessel's first or last segment coordinate is at the location of the node.
     */
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::RemoveAdjoiningVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{

    for(int i = 0; i < pAdjoiningVessels.size(); i++)
    {
        if (pAdjoiningVessels[i].lock() == vessel)
        {
            pAdjoiningVessels.erase(pAdjoiningVessels.begin() + i);
            i--;
        }
    }

    /*
     * todo perhaps should raise a warning if vessel is not attached to this node
     */

}


template<unsigned SPATIAL_DIM>
bool CaVascularNetworkNode<SPATIAL_DIM>::IsAttachedToVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{

    bool vesselIsAttachedToNode = false;

    for(int i = 0; i < pAdjoiningVessels.size(); i++)
    {
        if (pAdjoiningVessels[i].lock() == vessel)
        {
            assert(vessel->GetNode1() == shared() || vessel->GetNode2() == shared());
            vesselIsAttachedToNode = true;
            break;
        }
    }

    return vesselIsAttachedToNode;
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::SetIsInputNode(bool value)
{
    mIsInputNode = value;
}

template<unsigned SPATIAL_DIM>
void CaVascularNetworkNode<SPATIAL_DIM>::SetIsOutputNode(bool value)
{
    mIsOutputNode = value;
}

// Explicit instantiation

template class CaVascularNetworkNode<2>;
template class CaVascularNetworkNode<3>;

