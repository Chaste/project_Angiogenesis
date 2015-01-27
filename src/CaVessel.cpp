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

#include "CaVessel.hpp"

template<unsigned DIM>
CaVessel<DIM>::CaVessel(boost::shared_ptr<CaVesselSegment<DIM> > pSegment)
	: mSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >()),
	  mpDataContainer(boost::shared_ptr<VascularNetworkData>(new VascularNetworkData())),
	  mId(0),
	  mLabel("")
{
	mSegments.push_back(pSegment);
}

template<unsigned DIM>
CaVessel<DIM>::CaVessel(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments)
	: mSegments(segments),
	  mpDataContainer(boost::shared_ptr<VascularNetworkData>(new VascularNetworkData())),
	  mId(0),
	  mLabel("")
{
	for (unsigned i = 1; i < mSegments.size(); i++)
	{
		if(mSegments[i]->GetNodes(0) != mSegments[i-1]->GetNodes(1))
		{
			EXCEPTION("Input vessel segments are not attached in the correct order.");
		}
	}
}

template<unsigned DIM>
CaVessel<DIM>::~CaVessel()
{
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> >CaVessel<DIM>::Create(boost::shared_ptr<CaVesselSegment<DIM> > pSegment)
{
	boost::shared_ptr<CaVessel<DIM> > pSelf(new CaVessel<DIM>(pSegment));

	// Add the vessel to the segment
	pSegment->AddVessel(pSelf->shared_from_this());
	return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> >CaVessel<DIM>::Create(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments)
{
	boost::shared_ptr<CaVessel<DIM> > pSelf(new CaVessel<DIM>(segments));

	// Add the vessel to the segments
	for (unsigned i = 1; i < segments.size(); i++)
	{
		segments[i]->AddVessel(pSelf->shared_from_this());
	}
	return pSelf;
}

template<unsigned DIM>
void CaVessel<DIM>::AddSegments(boost::shared_ptr<CaVesselSegment<DIM> > pSegment)
{
	///\todo adding segments at the start could be slow. Maybe a deque is better than
	// vector in this case.

	if(pSegment->GetNodes(0) == mSegments.back()->GetNodes(1))
	// Append to end of vessel
	{
		mSegments.push_back(pSegment);
	}
	else if(pSegment->GetNodes(1) == mSegments.front()->GetNodes(0))
	// Insert at the start of the vessel
	{
		mSegments.insert(mSegments.begin(), pSegment);
	}
	else
	{
		EXCEPTION("Input vessel segment does not coincide with any end of the vessel.");
	}
}

template<unsigned DIM>
void CaVessel<DIM>::AddSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments)
{

	// Check that the first or end segment coincide with the start or end of the vessel
	bool start = false;
	bool end = false;

	if(segments.front()->GetNodes(0) == mSegments.back()->GetNodes(1))
	{
		start = true;
	}
	else if(segments.back()->GetNodes(1) == mSegments.front()->GetNodes(0))
	{
		end = true;
	}
	else
	{
		EXCEPTION("Input vessel segments do not coincide with any end of the vessel.");
	}

	// Check that the segments are in the correct order
	for (unsigned i = 1; i < segments.size(); i++)
	{
		if(segments[i]->GetNodes(0) != segments[i-1]->GetNodes(1))
		{
			EXCEPTION("Input vessel segments are not ordered correctly.");
		}
		segments[i]->AddVessel(Shared());
	}

	if (start)
	{
		mSegments.insert(mSegments.end(), segments.begin(), segments.end());
	}
	else if (end)
	{
		std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > temp = segments;
		temp.insert(temp.end(), mSegments.begin(), mSegments.end());
		mSegments = temp;
	}
}

template<unsigned DIM>
boost::shared_ptr<VascularNetworkData> CaVessel<DIM>::GetDataContainer()
{
	return mpDataContainer;
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetId()
{
	return mId;
}

template<unsigned DIM>
const std::string& CaVessel<DIM>::rGetLabel()
{
	return mLabel;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVessel<DIM>::GetStartNode()
{
    return mSegments.front()->GetNodes(0);
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVessel<DIM>::GetEndNode()
{
	return mSegments.back()->GetNodes(1);
}

template<unsigned DIM>
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > CaVessel<DIM>::GetSegments()
{
	return mSegments;
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVessel<DIM>::GetSegments(unsigned i)
{
	if(i < mSegments.size())
	{
		return mSegments[i];
	}
	else
	{
		EXCEPTION("Requested segment index exceeds number of segments.");
	}
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfSegments()
{
    return mSegments.size();
}

template<unsigned DIM>
void CaVessel<DIM>::SetDataContainer(boost::shared_ptr<VascularNetworkData> pDataContainer)
{
	mpDataContainer->SetMap(pDataContainer->GetMap());
}

template<unsigned DIM>
void CaVessel<DIM>::SetId(unsigned id)
{
	mId = id;
}

template<unsigned DIM>
void CaVessel<DIM>::SetLabel(const std::string& label)
{
	mLabel = label;
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVessel<DIM>::Shared()
{
    return this->shared_from_this();
}

// Explicit instantiation
template class CaVessel<2>;
template class CaVessel<3>;
