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

#include <algorithm>
#include <stddef.h>
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"

template<unsigned DIM>
VascularNode<DIM>::VascularNode(const ChastePoint<DIM>& rLocation)
	: mLocation(rLocation),
	  mpCell(CellPtr()),
	  mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
	  mDataContainer(),
	  mId(0),
	  mLabel(""),
	  mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
	  mPressure(0.0),
	  mRadius(0.0),
	  mIsInputNode(false),
	  mIsOutputNode(false),
	  mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(const VascularNode<DIM>& rExistingNode)
: boost::enable_shared_from_this<VascularNode<DIM> >(),
  mLocation(rExistingNode.GetLocation()),
  mpCell(CellPtr()),
  mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
  mDataContainer(),
  mId(0),
  mLabel(""),
  mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
  mPressure(rExistingNode.GetPressure()),
  mRadius(rExistingNode.GetRadius()),
  mIsInputNode(),
  mIsOutputNode(),
  mIsMigrating()
{
    mDataContainer.SetMap(rExistingNode.rGetDataContainer().GetMap());
    mIsInputNode = rExistingNode.IsInputNode();
    mIsOutputNode = rExistingNode.IsOutputNode();
    mIsMigrating = rExistingNode.IsMigrating();
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(double point1, double point2, double point3)
    : mLocation(ChastePoint<DIM>(point1, point2, point3)),
      mpCell(CellPtr()),
      mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
      mDataContainer(),
      mId(0),
      mLabel(""),
      mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
	  mPressure(0.0),
	  mRadius(1.0),
	  mIsInputNode(false),
	  mIsOutputNode(false),
      mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(c_vector<double, DIM> location)
    : mLocation(ChastePoint<DIM>(location)),
      mpCell(CellPtr()),
      mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
      mDataContainer(),
      mId(0),
      mLabel(""),
      mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
	  mPressure(0.0),
	  mRadius(1.0),
	  mIsInputNode(false),
	  mIsOutputNode(false),
      mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::~VascularNode()
{
}

template<unsigned DIM>
void VascularNode<DIM>::AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment)
{
	// Vessel segments can only be attached to a node once.
	int number_times_attached_to_node = 0;

	for(unsigned i = 0; i < mVesselSegments.size(); i++)
	{
		if (mVesselSegments[i].lock() == pVesselSegment)
		{
			number_times_attached_to_node++;
		}
	}

	if (number_times_attached_to_node == 0)
	{
		mVesselSegments.push_back(boost::weak_ptr<CaVesselSegment<DIM> > (pVesselSegment));
	}
	else
	{
		EXCEPTION("This segment is already attached to this node.");
	}
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> >VascularNode<DIM>::Create(const ChastePoint<DIM>& rLocation)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (rLocation));
	return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> >VascularNode<DIM>::Create(double point1, double point2, double point3)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (point1, point2, point3));
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> >VascularNode<DIM>::Create(c_vector<double, DIM> location)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (location));
    return pSelf;
}

template<unsigned DIM>
CellPtr VascularNode<DIM>::GetCell() const
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
const VasculatureData& VascularNode<DIM>::rGetDataContainer() const
{
	return mDataContainer;
}

template<unsigned DIM>
template<typename T> T VascularNode<DIM>::GetData(const std::string& variableName)
{
	return mDataContainer.GetData<T>(variableName);
}

template<unsigned DIM>
std::vector<std::string> VascularNode<DIM>::GetDataKeys(bool castable_to_double) const
{
	return mDataContainer.GetKeys(castable_to_double);
}

template<unsigned DIM>
double VascularNode<DIM>::GetDistance(const ChastePoint<DIM>& rPoint)
{
	double distance_squared;
	ChastePoint<DIM> location = GetLocation();

	distance_squared = pow(location[0] -rPoint[0],2) + pow(location[1] -rPoint[1],2);
	if(DIM==3)
	{
		distance_squared += pow(location[2] -rPoint[2],2);
	}

	return std::sqrt(distance_squared);
}

template<unsigned DIM>
unsigned VascularNode<DIM>::GetId() const
{
	return mId;
}

template<unsigned DIM>
const std::string& VascularNode<DIM>::rGetLabel() const
{
	return mLabel;
}

template<unsigned DIM>
ChastePoint<DIM> VascularNode<DIM>::GetLocation() const
{
	if(mpCell)
	{
		return mpCellPopulation->GetLocationOfCellCentre(mpCell);
	}
	else
	{
		return mLocation;
	}
}

template<unsigned DIM>
double VascularNode<DIM>::GetPressure() const
{
	return mPressure;
}

template<unsigned DIM>
double VascularNode<DIM>::GetRadius() const
{
	return mRadius;
}

template<unsigned DIM>
unsigned VascularNode<DIM>::GetNumberOfSegments() const
{
	return mVesselSegments.size();
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > VascularNode<DIM>::GetVesselSegment(unsigned index) const
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
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >VascularNode<DIM>::GetVesselSegments() const
{
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;

	for(unsigned i=0; i<mVesselSegments.size(); i++)
	{
		segments.push_back(mVesselSegments[i].lock());
	}

	return segments;
}

template<unsigned DIM>
bool VascularNode<DIM>::HasCell() const
{
	return mpCell;
}

template<unsigned DIM>
bool VascularNode<DIM>::HasDataKey(const std::string& rKey) const
{
	return mDataContainer.HasKey(rKey);
}

template<unsigned DIM>
bool VascularNode<DIM>::IsAttachedTo(const boost::shared_ptr<CaVesselSegment<DIM> > pSegment) const
{
	if (pSegment->GetNode(0) == this->shared_from_this() || pSegment->GetNode(1) == this->shared_from_this())
	{
		return true;
	}
	else
	{
		return false;
	}
}

template<unsigned DIM>
bool VascularNode<DIM>::IsCoincident(const ChastePoint<DIM>& rPoint) const
{
	return (GetLocation().IsSamePoint(rPoint));
}

template<unsigned DIM>
bool VascularNode<DIM>::IsCoincident(const boost::shared_ptr<VascularNode<DIM> > pNode) const
{
	return (GetLocation().IsSamePoint(pNode->GetLocation()));
}


template<unsigned DIM>
bool VascularNode<DIM>::IsInputNode() const
{
    return mIsInputNode;
}


template<unsigned DIM>
bool VascularNode<DIM>::IsOutputNode() const
{
    return mIsOutputNode;
}

template<unsigned DIM>
bool VascularNode<DIM>::IsMigrating() const
{
    return mIsMigrating;
}

template<unsigned DIM>
void VascularNode<DIM>::RemoveCell()
{
	mpCell = CellPtr();
}

template<unsigned DIM>
void VascularNode<DIM>::RemoveSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment)
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

template<unsigned DIM>
void VascularNode<DIM>::SetCell(CellPtr pCell)
{
	if(mpCellPopulation)
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
void VascularNode<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation)
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
void VascularNode<DIM>::SetDataContainer(const VasculatureData& rDataContainer)
{
	mDataContainer = rDataContainer;
}

template<unsigned DIM>
template<typename T> void VascularNode<DIM>::SetData(const std::string& variableName, T value)
{
	mDataContainer.SetData(variableName, value);
}

template<unsigned DIM>
void VascularNode<DIM>::SetId(unsigned id)
{
	mId = id;
}

template<unsigned DIM>
void VascularNode<DIM>::IsInputNode(bool inputNode)
{
    mIsInputNode = inputNode;
}


template<unsigned DIM>
void VascularNode<DIM>::IsOutputNode(bool outputNode)
{
    mIsOutputNode = outputNode;
}

template<unsigned DIM>
void VascularNode<DIM>::IsMigrating(bool isMigrating)
{
    mIsMigrating = isMigrating;
}

template<unsigned DIM>
void VascularNode<DIM>::SetLabel(const std::string& rLabel)
{
	mLabel = rLabel;
}

template<unsigned DIM>
void VascularNode<DIM>::SetLocation(const ChastePoint<DIM>& rLocation)
{
	if (mpCell)
	{
		mpCell = CellPtr();
	}
	mLocation = rLocation;
}

template<unsigned DIM>
void VascularNode<DIM>::SetPressure(double pressure)
{
	mPressure = pressure;
}

template<unsigned DIM>
void VascularNode<DIM>::SetRadius(double radius)
{
	mRadius = radius;
}

// Explicit instantiation
template class VascularNode<2>;
template class VascularNode<3>;

template bool VascularNode<2>::GetData<bool>(const std::string& variableName);
template double VascularNode<2>::GetData<double>(const std::string& variableName);
template unsigned VascularNode<2>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> VascularNode<2>::GetData<std::vector<double> >(const std::string& variableName);
template void VascularNode<2>::SetData(const std::string& variableName, bool value);
template void VascularNode<2>::SetData(const std::string& variableName, double value);
template void VascularNode<2>::SetData(const std::string& variableName, unsigned value);
template void VascularNode<2>::SetData(const std::string& variableName, std::vector<double> value);

template bool VascularNode<3>::GetData<bool>(const std::string& variableName);
template double VascularNode<3>::GetData<double>(const std::string& variableName);
template unsigned VascularNode<3>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> VascularNode<3>::GetData<std::vector<double> >(const std::string& variableName);
template void VascularNode<3>::SetData(const std::string& variableName, bool value);
template void VascularNode<3>::SetData(const std::string& variableName, double value);
template void VascularNode<3>::SetData(const std::string& variableName, unsigned value);
template void VascularNode<3>::SetData(const std::string& variableName, std::vector<double> value);


