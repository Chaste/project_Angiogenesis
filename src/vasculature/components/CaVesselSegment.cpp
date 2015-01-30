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

#include "CaVesselSegment.hpp"

template<unsigned DIM>
CaVesselSegment<DIM>::CaVesselSegment(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2)
	: mNodes(std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > (pNode1, pNode2)),
	  mpDataContainer(boost::shared_ptr<VasculatureData>(new VasculatureData())),
	  mId(0),
	  mLabel(""),
	  mVessel(boost::weak_ptr<CaVessel<DIM> >())
{
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> >CaVesselSegment<DIM>::Create(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2)
{
	boost::shared_ptr<CaVesselSegment<DIM> > pSelf(new CaVesselSegment<DIM>(pNode1, pNode2));

	// Make sure the specified nodes are not the same
	if (pNode1 == pNode2)
	{
		EXCEPTION("Attempted to assign the same node to both ends of a vessel segment.");
	}

	// Add the segment to the nodes
	pNode1->AddSegment(pSelf->shared_from_this());
	pNode2->AddSegment(pSelf->shared_from_this());
	return pSelf;
}

template<unsigned DIM>
CaVesselSegment<DIM>::~CaVesselSegment()
{
}

template<unsigned DIM>
boost::shared_ptr<VasculatureData> CaVesselSegment<DIM>::GetDataContainer()
{
	return mpDataContainer;
}

template<unsigned DIM>
unsigned CaVesselSegment<DIM>::GetId()
{
	return mId;
}

template<unsigned DIM>
const std::string& CaVesselSegment<DIM>::rGetLabel()
{
	return mLabel;
}
template<unsigned DIM>
double CaVesselSegment<DIM>::GetLength()
{
	ChastePoint<DIM> point1 = mNodes.first->GetLocation();
	ChastePoint<DIM> point2 = mNodes.second->GetLocation();

	double length_squared = pow(point2[0] - point1[0], 2) + pow(point2[1] - point1[1], 2);

	if (DIM == 3u)
	{
		length_squared += pow(point2[2] - point1[2], 2);
	}
	return std::sqrt(length_squared);
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVesselSegment<DIM>::GetVessel()
{
	if(mVessel.lock())
	{
		return mVessel.lock();
	}
	else
	{
		EXCEPTION("A vessel has been requested but this segment doesn't have one.");
	}
}

template<unsigned DIM>
std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > CaVesselSegment<DIM>::GetNodes()
{
	return mNodes;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVesselSegment<DIM>::GetNodes(unsigned index)
{
	if (index == 0u)
	{
		return mNodes.first;
	}
	else if (index == 1u)
	{
		return mNodes.second;
	}
	else
	{
		EXCEPTION("A node index other than 0 or 1 has been requested for a Vessel Segment.");
	}
}

template<unsigned DIM>
void CaVesselSegment<DIM>::ReplaceNode(unsigned old_node_index, boost::shared_ptr<VascularNode<DIM> >  pNewNode)
{
	if (old_node_index == 0u)
	{
		mNodes.first->RemoveSegment(Shared());
		mNodes.first = pNewNode;
		mNodes.first->AddSegment(Shared());
	}
	else if (old_node_index == 1u)
	{
		mNodes.second->RemoveSegment(Shared());
		mNodes.second = pNewNode;
		mNodes.second->AddSegment(Shared());
	}
	else
	{
		EXCEPTION("A node index other than 0 or 1 has been requested for a Vessel Segment.");
	}
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetDataContainer(boost::shared_ptr<VasculatureData> pDataContainer)
{
	mpDataContainer->SetMap(pDataContainer->GetMap());
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetId(unsigned id)
{
	mId = id;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetLabel(const std::string& label)
{
	mLabel = label;
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVesselSegment<DIM>::Shared()
{
		boost::shared_ptr<CaVesselSegment<DIM> > pSegment = this->shared_from_this();
		return pSegment;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel)
{
	if(mVessel.lock())
	{
		EXCEPTION("This segment already has a vessel.");
	}
	mVessel = pVessel;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::RemoveVessel()
{
	mVessel = boost::weak_ptr<CaVessel<DIM> >();
}

// Explicit instantiation
template class CaVesselSegment<2>;
template class CaVesselSegment<3>;
