/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "CaVascularNetworkNode.hpp"
#include "CaVessel.hpp"

template<unsigned DIM>
CaVascularNetworkNode<DIM>::CaVascularNetworkNode()
	: mLocation(),
	  mPressure(0.0),
	  mAdjoiningVessels(),
	  mIsInputNode(false),
	  mIsOutputNode(false)
{
}

template<unsigned DIM>
CaVascularNetworkNode<DIM>::~CaVascularNetworkNode()
{
}

///\todo This method will return a bad weak pointer exception if a shared pointer
// to the node has not been previously created. A generic catch has been used, it
// should be replaced with one specific to the pointer.
//
template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVascularNetworkNode<DIM>::shared()
{
	try
	{
		boost::shared_ptr<CaVascularNetworkNode<DIM> > pNode = this->shared_from_this();
		return pNode;
	}
	catch(...)
	{
		EXCEPTION("Failed to return a shared pointer to a node. The node should be owned by a shared pointer before calling shared.");
	}
}

template<unsigned DIM>
ChastePoint<DIM> CaVascularNetworkNode<DIM>::GetLocation()
{
    return mLocation;
}

template<unsigned DIM>
double CaVascularNetworkNode<DIM>::GetPressure()
{
    return mPressure;
}

template<unsigned DIM>
unsigned CaVascularNetworkNode<DIM>::GetNumberOfAdjoiningVessels()
{
    return mAdjoiningVessels.size();
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetworkNode<DIM>::GetAdjoiningVessel(unsigned i)
{
	if(i < mAdjoiningVessels.size())
	{
		return mAdjoiningVessels[i].lock(); // lock() converts weak pointer (stored here) to shared pointer
	}
	else
	{
		EXCEPTION("Attempted to access a vessel with an out of range index.");
	}
}

template<unsigned DIM>
bool CaVascularNetworkNode<DIM>::IsInputNode()
{
    return mIsInputNode;
}

template<unsigned DIM>
bool CaVascularNetworkNode<DIM>::IsOutputNode()
{
    return mIsOutputNode;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetLocation(ChastePoint<DIM> location)
{
    mLocation = location;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetPressure(double pressure)
{
    mPressure = pressure;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::AddAdjoiningVessel(boost::shared_ptr<CaVessel<DIM> > vessel)
{
    // The same vessel may be adjoint to a single node twice but that node must be both node1 and node2 of the doubly adjoint vessel
    int number_times_attached_to_node = 0;

    for(unsigned i = 0; i < mAdjoiningVessels.size(); i++)
    {
        if (mAdjoiningVessels[i].lock() == vessel)
        {
        	number_times_attached_to_node++;
        }
    }

    if (number_times_attached_to_node == 0)
    {
        mAdjoiningVessels.push_back(boost::weak_ptr<CaVessel<DIM> >(vessel));
    }
    else if (number_times_attached_to_node == 1 && vessel->GetNode1() == shared() && vessel->GetNode2() == shared())
    {
        mAdjoiningVessels.push_back(boost::weak_ptr<CaVessel<DIM> >(vessel));
    }
    else
    {
        if (number_times_attached_to_node == 2)
        {
            EXCEPTION("Vessel is already attached to node twice (at both ends). Cannot attach vessel to same node again.");
        }
        else
        {
        	EXCEPTION("Vessels and nodes in inconsistent state.");
        }
    }

    /*
        todo should check that vessel being adjoined to node has not got an actively migrating tip located at this node
        todo should perhaps check that the vessel exists at the location where the node is located
        - i.e. that the vessel's first or last segment coordinate is at the location of the node.
     */
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::RemoveAdjoiningVessel(boost::shared_ptr<CaVessel<DIM> > vessel)
{
    bool attached_to_vessel = false;
    for(unsigned i = 0; i < mAdjoiningVessels.size(); i++)
    {
        if (mAdjoiningVessels[i].lock() == vessel)
        {
        	attached_to_vessel = true;
            mAdjoiningVessels.erase(mAdjoiningVessels.begin() + i);
            i--;
        }
    }

    if (!attached_to_vessel)
    {
    	EXCEPTION("Attempted to remove a vessel from a node it is not attached to.");
    }
}


template<unsigned DIM>
bool CaVascularNetworkNode<DIM>::IsAttachedToVessel(boost::shared_ptr<CaVessel<DIM> > vessel)
{
    bool attached_to_node = false;

    for(unsigned i = 0; i < mAdjoiningVessels.size(); i++)
    {
        if (mAdjoiningVessels[i].lock() == vessel)
        {
            if (vessel->GetNode1() != shared() && vessel->GetNode2() != shared())
            {
            	EXCEPTION("The vessel has been added to the node, but the node has not been added to the vessel.");
            }
            attached_to_node = true;
            break;
        }
    }
    return attached_to_node;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetIsInputNode(bool value)
{
    mIsInputNode = value;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetIsOutputNode(bool value)
{
    mIsOutputNode = value;
}

// Explicit instantiation
template class CaVascularNetworkNode<2>;
template class CaVascularNetworkNode<3>;

