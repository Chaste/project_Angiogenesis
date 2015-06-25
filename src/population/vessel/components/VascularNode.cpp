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

#include <boost/weak_ptr.hpp>
#include "UblasIncludes.hpp"
#include "Exception.hpp"
#include "SmartPointers.hpp"
#include "CaVesselSegment.hpp"

#include "VascularNode.hpp"

template<unsigned DIM>
VascularNode<DIM>::VascularNode(const ChastePoint<DIM>& rLocation) :
        mLocation(rLocation),
        mpCell(CellPtr()),
        mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
        mDataContainer(),
        mId(0),
        mLabel(),
        mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
        mRadius(0.0),
        mpNodeFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(double v1, double v2, double v3) :
        mLocation(ChastePoint<DIM>(v1, v2, v3)),
        mpCell(CellPtr()),
        mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
        mDataContainer(),
        mId(0),
        mLabel(),
        mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
        mRadius(1.0),
        mpNodeFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(c_vector<double, DIM> location) :
        mLocation(ChastePoint<DIM>(location)),
        mpCell(CellPtr()),
        mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
        mDataContainer(),
        mId(0),
        mLabel(),
        mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
        mRadius(1.0),
        mpNodeFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(const VascularNode<DIM>& rExistingNode) :
        boost::enable_shared_from_this<VascularNode<DIM> >(),
        mLocation(rExistingNode.GetLocation()),
        mpCell(CellPtr()),
        mpCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> >()),
        mDataContainer(),
        mId(0),
        mLabel(),
        mVesselSegments(std::vector<boost::weak_ptr<CaVesselSegment<DIM> > >()),
        mRadius(rExistingNode.GetRadius()),
        mpNodeFlowProperties(boost::shared_ptr<NodeFlowProperties>()),
        mIsMigrating()
{
    mDataContainer.SetMap(rExistingNode.rGetDataContainer().GetMap());
    SetFlowProperties(*(rExistingNode.GetFlowProperties()));
    mIsMigrating = rExistingNode.IsMigrating();
}


template<unsigned DIM>
VascularNode<DIM>::~VascularNode()
{
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNode<DIM>::Create(const ChastePoint<DIM>& rLocation)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (rLocation));
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNode<DIM>::Create(double v1, double v2, double v3)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (v1, v2, v3));
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNode<DIM>::Create(c_vector<double, DIM> location)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (location));
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNode<DIM>::Create(const VascularNode<DIM>& rExistingNode)
{
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (rExistingNode));
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNode<DIM>::Create(boost::shared_ptr<VascularNode<DIM> > pExistingNode)
{
    if(!pExistingNode)
    {
        EXCEPTION("A Null pointer cannot be used when copying nodes.");
    }
    MAKE_PTR_ARGS(VascularNode<DIM>, pSelf, (*pExistingNode));
    return pSelf;
}

template<unsigned DIM>
void VascularNode<DIM>::AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment)
{
    // Vessel segments can only be attached to a node once.
    for (unsigned idx = 0; idx < mVesselSegments.size(); idx++)
    {
        if (mVesselSegments[idx].lock() == pVesselSegment)
        {
            EXCEPTION("This segment is already attached to this node.");
        }
    }
    mVesselSegments.push_back(boost::weak_ptr<CaVesselSegment<DIM> >(pVesselSegment));
}

template<unsigned DIM>
CellPtr VascularNode<DIM>::GetCell() const
{
    if (!mpCell)
    {
        EXCEPTION("A Cell has been requested but none have been assigned to this Node.");
    }
    return mpCell;
}

template<unsigned DIM>
const VasculatureData& VascularNode<DIM>::rGetDataContainer() const
{
    return mDataContainer;
}

template<unsigned DIM>
template<typename T> T VascularNode<DIM>::GetData(const std::string& rKey) const
{
    return mDataContainer.GetData<T>(rKey);
}

template<unsigned DIM>
std::vector<std::string> VascularNode<DIM>::GetDataKeys(bool castableToDouble) const
{
    return mDataContainer.GetKeys(castableToDouble);
}

template<unsigned DIM>
double VascularNode<DIM>::GetDistance(boost::shared_ptr<VascularNode<DIM> > pNode) const
{
    return norm_2(pNode->GetLocation().rGetLocation() - mLocation.rGetLocation());
}

template<unsigned DIM>
double VascularNode<DIM>::GetDistance(const ChastePoint<DIM>& rPoint) const
{
    return norm_2(rPoint.rGetLocation() - mLocation.rGetLocation());
}

template<unsigned DIM>
double VascularNode<DIM>::GetDistance(c_vector<double, DIM> location) const
{
    return norm_2(location - mLocation.rGetLocation());
}

template<unsigned DIM>
boost::shared_ptr<NodeFlowProperties> VascularNode<DIM>::GetFlowProperties() const
{
    return mpNodeFlowProperties;
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
//    if (mpCell)
//    {
//        return mpCellPopulation->GetLocationOfCellCentre(mpCell);
//    }
//    else
//    {
        return mLocation;
//    }
}

template<unsigned DIM>
c_vector<double, DIM> VascularNode<DIM>::GetLocationVector() const
{
//    if (mpCell)
//    {
//        return mpCellPopulation->GetLocationOfCellCentre(mpCell);
//    }
//    else
//    {
        return mLocation.rGetLocation();
//    }
}

template<unsigned DIM>
unsigned VascularNode<DIM>::GetNumberOfSegments() const
{
    return mVesselSegments.size();
}

template<unsigned DIM>
double VascularNode<DIM>::GetRadius() const
{
    return mRadius;
}

template<unsigned DIM>
std::map<std::string, double> VascularNode<DIM>::GetVtkData() const
{
    std::map<std::string, double> vtk_data;
    std::map<std::string, double> flow_data = mpNodeFlowProperties->GetVtkData();
    vtk_data.insert(flow_data.begin(), flow_data.end());
    vtk_data["Node Id"] = double(GetId());
    vtk_data["Node Radius"] = GetRadius();
    vtk_data["Node Is Migrating"] = double(IsMigrating());
    return vtk_data;
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > VascularNode<DIM>::GetVesselSegment(unsigned index) const
{
    if (index >= mVesselSegments.size())
    {
        EXCEPTION("Attempted to access a segment with an out of range index.");
    }
    return mVesselSegments[index].lock();
}

template<unsigned DIM>
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > VascularNode<DIM>::GetVesselSegments() const
{
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments(mVesselSegments.size());
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        segments[idx] = mVesselSegments[idx].lock();
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
    bool is_attached = false;
    if (pSegment->GetNode(0) == this->shared_from_this() || pSegment->GetNode(1) == this->shared_from_this())
    {
        is_attached = true;
    }
    return is_attached;
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
    for (unsigned idx = 0; idx < mVesselSegments.size(); idx++)
    {
        if (mVesselSegments[idx].lock() == pVesselSegment)
        {
            mVesselSegments.erase(mVesselSegments.begin() + idx);
            break;
        }
    }
}

template<unsigned DIM>
void VascularNode<DIM>::SetCell(CellPtr pCell)
{
//    if (!mpCellPopulation)
//    {
//        EXCEPTION("Attempted to add a Cell without first adding a CellPopulation.");
//    }
//
//    std::list<CellPtr> cell_list = mpCellPopulation->rGetCells();
//    bool found = (std::find(cell_list.begin(), cell_list.end(), pCell) != cell_list.end());
//    if (!found)
//    {
//        EXCEPTION("Attempted to add a Cell that is not in the assigned CellPopulation.");
//    }

    mpCell = pCell;
}

template<unsigned DIM>
void VascularNode<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation)
{
    // If there is an existing cell not in the cell population then remove it.
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
template<typename T> void VascularNode<DIM>::SetData(const std::string& rKey, T value)
{
    mDataContainer.SetData(rKey, value);
}

template<unsigned DIM>
void VascularNode<DIM>::SetFlowProperties(const NodeFlowProperties& rFlowProperties)
{
    mpNodeFlowProperties = boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties(rFlowProperties));
}

template<unsigned DIM>
void VascularNode<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void VascularNode<DIM>::SetIsMigrating(bool isMigrating)
{
    if (isMigrating)
    {
        if (GetNumberOfSegments() > 1)
        {
            EXCEPTION("A node should not migrate if its connectivity > 1");
        }
    }
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
void VascularNode<DIM>::SetLocation(c_vector<double, DIM> location)
{
    if (mpCell)
    {
        mpCell = CellPtr();
    }
    mLocation = ChastePoint<DIM>(location);
}

template<unsigned DIM>
void VascularNode<DIM>::SetRadius(double radius)
{
    mRadius = radius;
}

// Explicit instantiation
template class VascularNode<2> ;
template class VascularNode<3> ;

template bool VascularNode<2>::GetData<bool>(const std::string& rKey) const;
template double VascularNode<2>::GetData<double>(const std::string& rKey) const;
template unsigned VascularNode<2>::GetData<unsigned>(const std::string& rKey) const;
template std::vector<double> VascularNode<2>::GetData<std::vector<double> >(const std::string& rKey) const;
template void VascularNode<2>::SetData(const std::string& rKey, bool value);
template void VascularNode<2>::SetData(const std::string& rKey, double value);
template void VascularNode<2>::SetData(const std::string& rKey, unsigned value);
template void VascularNode<2>::SetData(const std::string& rKey, std::vector<double> value);

template bool VascularNode<3>::GetData<bool>(const std::string& rKey) const;
template double VascularNode<3>::GetData<double>(const std::string& rKey) const;
template unsigned VascularNode<3>::GetData<unsigned>(const std::string& rKey) const;
template std::vector<double> VascularNode<3>::GetData<std::vector<double> >(const std::string& rKey) const;
template void VascularNode<3>::SetData(const std::string& rKey, bool value);
template void VascularNode<3>::SetData(const std::string& rKey, double value);
template void VascularNode<3>::SetData(const std::string& rKey, unsigned value);
template void VascularNode<3>::SetData(const std::string& rKey, std::vector<double> value);
