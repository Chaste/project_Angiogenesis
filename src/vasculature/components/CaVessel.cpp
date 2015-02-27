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
  mDataContainer(),
  mId(0),
  mLabel("")
  {
	mSegments.push_back(pSegment);
  }


template<unsigned DIM>
CaVessel<DIM>::CaVessel(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments)
: mSegments(segments),
  mDataContainer(),
  mId(0),
  mLabel("")
  {

	if (segments.size() > 1)
	{
        for (unsigned i = 1; i < mSegments.size(); i++)
        {
            if(!mSegments[i]->IsConnectedTo(mSegments[i-1]))
            {
                EXCEPTION("Input vessel segments are not attached in the correct order.");
            }
        }

        for (unsigned i = 0; i < mSegments.size(); i++)
        {
            for (unsigned j = 0; j < mSegments.size(); j++)
            {
                if (i != j && i != j-1 && i != j+1)
                {
                    if(mSegments[i]->IsConnectedTo(mSegments[j]))
                    {
                        EXCEPTION("Input vessel segments are not correctly connected.");
                    }
                }
            }
        }
    }
  }
  
template<unsigned DIM>
CaVessel<DIM>::CaVessel(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes)
: mSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >()),
  mDataContainer(),
  mId(0),
  mLabel("")
  {

	if (nodes.size() < 2)
	{
		EXCEPTION("Insufficient number of nodes to define a segment.");
	}
	else
	{
        for (unsigned i = 1; i < nodes.size(); i++)
        {
        	MAKE_VN_PTR_ARGS(CaVesselSegment<DIM>, pSegment, (nodes[i-1], nodes[i]));
        	mSegments.push_back(pSegment);
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
	for (unsigned i = 0; i < segments.size(); i++)
	{
		segments[i]->AddVessel(pSelf->shared_from_this());
	}
	return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> >CaVessel<DIM>::Create(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes)
{
	boost::shared_ptr<CaVessel<DIM> > pSelf(new CaVessel<DIM>(nodes));

	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = pSelf->GetSegments();

	// Add the vessel to the new segments
	for (unsigned i = 0; i < segments.size(); i++)
	{
		segments[i]->AddVessel(pSelf->shared_from_this());
	}
	return pSelf;
}

template<unsigned DIM>
void CaVessel<DIM>::AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pSegment)
{

	if(pSegment->IsConnectedTo(mSegments.back()))
		// Append to end of vessel
	{
		pSegment->AddVessel(Shared());
		mSegments.push_back(pSegment);
	}
	else if(pSegment->IsConnectedTo(mSegments.front()))
		// Insert at the start of the vessel
	{
		pSegment->AddVessel(Shared());
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

	if(segments.front()->IsConnectedTo(mSegments.back()))
	{
		mSegments.insert(mSegments.end(), segments.begin(), segments.end());
	}
	else if(segments.back()->IsConnectedTo(mSegments.front()))
	{
		segments.insert(segments.end(), mSegments.begin(), mSegments.end());
		mSegments = segments;
	}
	else if(segments.front()->IsConnectedTo(mSegments.front()))
	{
		std::reverse(segments.begin(),segments.end());
		segments.insert(segments.end(), mSegments.begin(), mSegments.end());
		mSegments = segments;
	}
	else if(segments.back()->IsConnectedTo(mSegments.back()))
	{
		std::reverse(segments.begin(),segments.end());
		mSegments.insert(mSegments.end(), segments.begin(), segments.end());
	}
	else
	{
		EXCEPTION("Input vessel segments do not coincide with any end of the vessel.");
	}

	for (unsigned i = 1; i < mSegments.size(); i++)
	{
		if(!mSegments[i]->IsConnectedTo(mSegments[i-1]))
		{
			EXCEPTION("Input vessel segments are not attached in the correct order.");
		}
	}

	for (unsigned i = 0; i < mSegments.size(); i++)
	{
		for (unsigned j = 0; j < mSegments.size(); j++)
		{
			if (i != j && i != j-1 && i != j+1)
			{
				if(mSegments[i]->IsConnectedTo(mSegments[j]))
				{
					EXCEPTION("Input vessel segments are not correctly connected.");
				}
			}
		}
	}

	// Add the vessel to the segments
	for (unsigned i = 0; i < segments.size(); i++)
	{
		segments[i]->AddVessel(Shared());
	}

}

template<unsigned DIM>
const VasculatureData&  CaVessel<DIM>::rGetDataContainer() const
{
	return mDataContainer;
}

template<unsigned DIM>
std::vector<std::string> CaVessel<DIM>::GetDataKeys(bool castable_to_double) const
{
	return mDataContainer.GetKeys(castable_to_double);
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVessel<DIM>::GetEndNode()
{
	if (mSegments.size() == 1)
	{
		return mSegments.back()->GetNode(1);
	}
	else
	{
		if (mSegments[mSegments.size() - 2]->HasNode(mSegments.back()->GetNode(0)))
		{
			return mSegments.back()->GetNode(1);
		}
		else if (mSegments[mSegments.size() - 2]->HasNode(mSegments.back()->GetNode(1)))
		{
			return mSegments.back()->GetNode(0);
		}
		else
		{
			EXCEPTION("Vessel segments at end of vessel are not connected correctly.");
		}
	}
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
double CaVessel<DIM>::GetLength()
{
	double length = 0.0;

	for (unsigned i = 0; i < mSegments.size(); i++)
	{
		length += mSegments[i]->GetLength();
	}

	return length;
}

template <unsigned DIM>
std::set<boost::shared_ptr<VascularNode<DIM> > > CaVessel<DIM>::GetNodes()
{
	std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator jt;
	for(jt = mSegments.begin(); jt != mSegments.end(); jt++)
	{
		nodes.insert((*jt)->GetNode(0));
		nodes.insert((*jt)->GetNode(1));
	}

	return nodes;
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfNodes()
{
	return GetNumberOfSegments() + 1;
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfSegments()
{
	return mSegments.size();
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVessel<DIM>::GetSegment(unsigned i)
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
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > CaVessel<DIM>::GetSegments()
{
	return mSegments;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVessel<DIM>::GetStartNode()
{

	if (mSegments.size() == 1)
	{
		return mSegments.front()->GetNode(0);
	}
	else
	{
		if (mSegments[1]->HasNode(mSegments.front()->GetNode(0)))
		{
			return mSegments.front()->GetNode(1);
		}
		else if (mSegments[1]->HasNode(mSegments.front()->GetNode(1)))
		{
			return mSegments.front()->GetNode(0);
		}
		else
		{
			EXCEPTION("Vessel segments at start of vessel are not connected correctly.");
		}
	}
}

template<unsigned DIM>
bool CaVessel<DIM>::HasDataKey(const std::string& rKey) const
{
	return mDataContainer.HasKey(rKey);
}

template<unsigned DIM>
bool CaVessel<DIM>::IsConnectedTo(boost::shared_ptr<CaVessel<DIM> > pOtherVessel)
{
	if(GetStartNode() == pOtherVessel->GetStartNode() || GetEndNode() == pOtherVessel->GetStartNode() ||
			GetStartNode() == pOtherVessel->GetEndNode() || GetEndNode() == pOtherVessel->GetEndNode())
	{
		return true;
	}
	else
	{
		return false;
	}
}

template<unsigned DIM>
void CaVessel<DIM>::RemoveSegments(SegmentLocation::Value location)
{

	if(mSegments.size() == 1)
	{
		EXCEPTION("Vessel must have at least one segment.");
	}
	if (location == SegmentLocation::Start)
	{
		mSegments.front()->RemoveVessel();
		mSegments.erase(mSegments.begin());
	}
	else if (location == SegmentLocation::End)
	{
		mSegments.back()->RemoveVessel();
		mSegments.pop_back();
	}
	else
	{
		EXCEPTION("You can only remove segments from the start or end of vessels.");
	}
}

template<unsigned DIM>
void CaVessel<DIM>::SetDataContainer(const VasculatureData& rDataContainer)
{
	mDataContainer = rDataContainer;
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
