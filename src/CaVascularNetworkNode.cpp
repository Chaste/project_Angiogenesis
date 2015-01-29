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

template<unsigned DIM>
CaVascularNetworkNode<DIM>::CaVascularNetworkNode(ChastePoint<DIM> location)
	: mLocation(location),
	  mpCell(CellPtr()),
	  mpCellPopulation(NULL),
	  mpDataContainer(boost::shared_ptr<VascularNetworkData>(new VascularNetworkData())),
	  mId(0),
	  mLabel(""),
	  mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >())
{
}

template<unsigned DIM>
CaVascularNetworkNode<DIM>::~CaVascularNetworkNode()
{
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> >CaVascularNetworkNode<DIM>::Create(ChastePoint<DIM> location)
{
	boost::shared_ptr<CaVascularNetworkNode<DIM> > pSelf(new CaVascularNetworkNode<DIM>(location));
	return pSelf;
}

template<unsigned DIM>
CellPtr CaVascularNetworkNode<DIM>::GetCell()
{
	if(mpCell)
	{
		return mpCell;
	}
	else
	{
		EXCEPTION("A Cell has been requested but none have been assigned to this Node.");
	}
}

template<unsigned DIM>
unsigned CaVascularNetworkNode<DIM>::GetId()
{
	return mId;
}

template<unsigned DIM>
const std::string& CaVascularNetworkNode<DIM>::rGetLabel()
{
	return mLabel;
}

template<unsigned DIM>
ChastePoint<DIM> CaVascularNetworkNode<DIM>::GetLocation()
{
	if(mpCell)
	{
		ChastePoint<DIM> location(mpCellPopulation->GetLocationOfCellCentre(mpCell));
		return location;
	}
	else
	{
		return mLocation;
	}
}

template<unsigned DIM>
boost::shared_ptr<VascularNetworkData> CaVascularNetworkNode<DIM>::GetDataContainer()
{
	return mpDataContainer;
}

template<unsigned DIM>
bool CaVascularNetworkNode<DIM>::HasCell()
{
	return mpCell;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetDataContainer(boost::shared_ptr<VascularNetworkData> pDataContainer)
{
	mpDataContainer->SetMap(pDataContainer->GetMap());
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetCell(CellPtr pCell)
{
	if(mpCellPopulation != NULL)
	{
		std::list<CellPtr> cell_list = mpCellPopulation->rGetCells();
		bool found = (std::find(cell_list.begin(), cell_list.end(), pCell) != cell_list.end());
		if (found)
		{
			mpCell = pCell;
		}
		else
		{
			EXCEPTION("Attempted to add a Cell that is not in the assigned CellPopulation.");
		}
	}
	else
	{
		EXCEPTION("Attempted to add a Cell without first adding a CellPopulation.");
	}
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetCellPopulation(CaBasedCellPopulation<DIM>* pCellPopulation)
{
	if (mpCell)
	{
		std::list<CellPtr> cell_list = pCellPopulation->rGetCells();
		bool found = (std::find(cell_list.begin(), cell_list.end(), mpCell) != cell_list.end());

		if (!found)
		{
			mpCell = CellPtr();
		}
	}
	mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetId(unsigned id)
{
	mId = id;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetLabel(const std::string& label)
{
	mLabel = label;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::SetLocation(ChastePoint<DIM> location)
{
	if (mpCell)
	{
		mpCell = CellPtr();
	}
	mLocation = location;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::RemoveCell()
{
	mpCell = CellPtr();
}

template<unsigned DIM>
unsigned CaVascularNetworkNode<DIM>::GetNumberOfSegments()
{
	return mVesselSegments.size();
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVascularNetworkNode<DIM>::GetVesselSegments(unsigned index)
{
	if(index < mVesselSegments.size())
	{
		return mVesselSegments[index].lock();
	}
	else
	{
		EXCEPTION("Attempted to access a segment with an out of range index.");
	}
}

template<unsigned DIM>
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >CaVascularNetworkNode<DIM>::GetVesselSegments()
{
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;

	for(unsigned i=0; i<mVesselSegments.size(); i++)
	{
		segments.push_back(mVesselSegments[i].lock());
	}

	return segments;
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > vessel_segment)
{
	// Vessel segments can only be attached to a node once.
	int number_times_attached_to_node = 0;

	for(unsigned i = 0; i < mVesselSegments.size(); i++)
	{
		if (mVesselSegments[i].lock() == vessel_segment)
		{
			number_times_attached_to_node++;
		}
	}

	if (number_times_attached_to_node == 0)
	{
		mVesselSegments.push_back(boost::weak_ptr<CaVesselSegment<DIM> > (vessel_segment));
	}
	else
	{
		EXCEPTION("This segment is already attached to this node.");
	}
}

template<unsigned DIM>
void CaVascularNetworkNode<DIM>::RemoveSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment)
{
	for(unsigned i = 0; i < mVesselSegments.size(); i++)
	{
		if (mVesselSegments[i].lock() == pVesselSegment)
		{
			mVesselSegments.erase(mVesselSegments.begin() + i);
			i--;
		}
	}
}

// Explicit instantiation
template class CaVascularNetworkNode<2>;
template class CaVascularNetworkNode<3>;

