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

#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VascularNode.hpp"
#include "CaVessel.hpp"

#include "CaVesselSegment.hpp"

template<unsigned DIM>
CaVesselSegment<DIM>::CaVesselSegment(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2) :
        mNodes(std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > >(pNode1, pNode2)),
        mDataContainer(),
        mId(0),
        mLabel(),
        mVessel(boost::weak_ptr<CaVessel<DIM> >()),
        mRadius(1.0),
        mpFlowProperties(boost::shared_ptr<SegmentFlowProperties> (new SegmentFlowProperties()))
{
}

template<unsigned DIM>
CaVesselSegment<DIM>::CaVesselSegment(const CaVesselSegment<DIM>& rSegment) :
    boost::enable_shared_from_this<CaVesselSegment<DIM> >(),
    mNodes(rSegment.GetNodes()),
    mDataContainer(),
    mId(0),
    mLabel(),
    mVessel(boost::weak_ptr<CaVessel<DIM> >()),
    mRadius(rSegment.GetRadius()),
    mpFlowProperties(boost::shared_ptr<SegmentFlowProperties> ())
{
    mDataContainer.SetMap(rSegment.rGetDataContainer().GetMap());
    SetFlowProperties(*(rSegment.GetFlowProperties()));
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVesselSegment<DIM>::Create(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2)
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
boost::shared_ptr<CaVesselSegment<DIM> > CaVesselSegment<DIM>::Create(boost::shared_ptr<CaVesselSegment<DIM> > pSegment)
{
    if(!pSegment)
    {
        EXCEPTION("A Null pointer cannot be used when copying segments.");
    }
    MAKE_PTR_ARGS(CaVesselSegment<DIM>, pSelf, (*pSegment));

    // Add the segment to the nodes
    pSelf->GetNode(0)->AddSegment(pSelf->shared_from_this());
    pSelf->GetNode(1)->AddSegment(pSelf->shared_from_this());
    return pSelf;
}

template<unsigned DIM>
CaVesselSegment<DIM>::~CaVesselSegment()
{
}

template<unsigned DIM>
void CaVesselSegment<DIM>::AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel)
{
    mVessel = pVessel;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::CopyDataFromExistingSegment(const boost::shared_ptr<CaVesselSegment<DIM> > pTargetSegment)
{
    mDataContainer.SetMap(pTargetSegment->rGetDataContainer().GetMap());
    mRadius = pTargetSegment->GetRadius();
    SetFlowProperties(*(pTargetSegment->GetFlowProperties()));
}

template<unsigned DIM>
const VasculatureData& CaVesselSegment<DIM>::rGetDataContainer() const
{
    return mDataContainer;
}

template<unsigned DIM>
template<typename T> T CaVesselSegment<DIM>::GetData(const std::string& rKey) const
{
    return mDataContainer.GetData<T>(rKey);
}

template<unsigned DIM>
std::vector<std::string> CaVesselSegment<DIM>::GetDataKeys(bool castableToDouble) const
{
    return mDataContainer.GetKeys(castableToDouble);
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetDistance(c_vector<double, DIM> location) const
{
    c_vector<double, DIM> start_location = GetNode(0)->GetLocation().rGetLocation();
    c_vector<double, DIM> end_location = GetNode(1)->GetLocation().rGetLocation();

    c_vector<double, DIM> segment_vector = end_location - start_location;
    c_vector<double, DIM> point_vector = location - start_location;

    double dp_segment_point = inner_prod(segment_vector, point_vector);
    // Point projection is outside segment, return node0 distance
    if (dp_segment_point <= 0.0)
    {
        return GetNode(0)->GetDistance(location);
    }

    double dp_segment_segment = inner_prod(segment_vector, segment_vector);
    // Point projection is outside segment, return node1 distance
    if (dp_segment_segment <= dp_segment_point)
    {
        return GetNode(1)->GetDistance(location);
    }

    // Point projection is inside segment, get distance to point projection
    double projection_ratio = dp_segment_point / dp_segment_segment;
    c_vector<double, DIM> projected_point = start_location + projection_ratio * segment_vector;
    double distance = norm_2(projected_point - location);

    return distance;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetDistance(const ChastePoint<DIM>& rPoint) const
{
    return GetDistance(rPoint.rGetLocation());
}
template<unsigned DIM>
boost::shared_ptr<SegmentFlowProperties> CaVesselSegment<DIM>::GetFlowProperties() const
{
    return mpFlowProperties;
}

template<unsigned DIM>
unsigned CaVesselSegment<DIM>::GetId() const
{
    return mId;
}

template<unsigned DIM>
const std::string& CaVesselSegment<DIM>::rGetLabel() const
{
    return mLabel;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetLength() const
{
    return norm_2(mNodes.second->GetLocationVector() - mNodes.first->GetLocationVector());
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetRadius() const
{
    return mRadius;
}

template<unsigned DIM>
c_vector<double, DIM> CaVesselSegment<DIM>::GetMidPoint() const
{
    return (mNodes.second->GetLocationVector() + mNodes.first->GetLocationVector()) / 2.0;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVesselSegment<DIM>::GetNode(unsigned index) const
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
boost::shared_ptr<VascularNode<DIM> > CaVesselSegment<DIM>::GetOppositeNode(boost::shared_ptr<VascularNode<DIM> > pInputNode) const
{
    if(pInputNode == mNodes.first)
    {
        return mNodes.second;
    }
    else if(pInputNode == mNodes.second)
    {
        return mNodes.first;
    }
    else
    {
        EXCEPTION("Input node is not on the segment");
    }
}

template<unsigned DIM>
std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > CaVesselSegment<DIM>::GetNodes() const
{
    return mNodes;
}

template<unsigned DIM>
c_vector<double, DIM> CaVesselSegment<DIM>::GetPointProjection(c_vector<double, DIM> location) const
{
    c_vector<double, DIM> start_location = GetNode(0)->GetLocation().rGetLocation();
    c_vector<double, DIM> end_location = GetNode(1)->GetLocation().rGetLocation();

    c_vector<double, DIM> segment_vector = end_location - start_location;
    c_vector<double, DIM> point_vector = location - start_location;

    double dp_segment_point = inner_prod(segment_vector, point_vector);
    double dp_segment_segment = inner_prod(segment_vector, segment_vector);

    if (dp_segment_point <= 0.0 || dp_segment_segment <= dp_segment_point)
    {
        EXCEPTION("Projection of point is outside segment.");
    }

    // Point projection is inside segment, get distance to point projection
    double projection_ratio = dp_segment_point / dp_segment_segment;
    c_vector<double, DIM> projected_point = start_location + projection_ratio * segment_vector;
    return projected_point;
}

template<unsigned DIM>
c_vector<double, DIM> CaVesselSegment<DIM>::GetPointProjection(const ChastePoint<DIM>& rPoint) const
{
    return GetPointProjection(rPoint.rGetLocation());
}

template<unsigned DIM>
c_vector<double, DIM> CaVesselSegment<DIM>::GetUnitTangent() const
{
    return (mNodes.second->GetLocationVector() - mNodes.first->GetLocationVector()) / GetLength();
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVesselSegment<DIM>::GetVessel() const
{
    if (mVessel.lock())
    {
        return mVessel.lock();
    }
    else
    {
        EXCEPTION("A vessel has been requested but this segment doesn't have one.");
    }
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::HasDataKey(const std::string& rKey) const
{
    return mDataContainer.HasKey(rKey);
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::HasNode(boost::shared_ptr<VascularNode<DIM> > pNode) const
{
    return (pNode == GetNode(0) || pNode == GetNode(1));
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::IsConnectedTo(boost::shared_ptr<CaVesselSegment<DIM> > otherSegment) const
{
    bool isConnectedToSegment = false;
    if (this->GetNode(0) == otherSegment->GetNode(0) || this->GetNode(0) == otherSegment->GetNode(1)
            || this->GetNode(1) == otherSegment->GetNode(0) || this->GetNode(1) == otherSegment->GetNode(1))
    {
        isConnectedToSegment = true;
    }

    return isConnectedToSegment;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::RemoveVessel()
{
    mVessel = boost::weak_ptr<CaVessel<DIM> >();
}

template<unsigned DIM>
void CaVesselSegment<DIM>::ReplaceNode(unsigned oldNodeIndex, boost::shared_ptr<VascularNode<DIM> > pNewNode)
{
    if (oldNodeIndex == 0u)
    {
        mNodes.first->RemoveSegment(Shared());
        mNodes.first = pNewNode;
        mNodes.first->AddSegment(Shared());
    }
    else if (oldNodeIndex == 1u)
    {
        mNodes.second->RemoveSegment(Shared());
        mNodes.second = pNewNode;
        mNodes.second->AddSegment(Shared());
    }
    else
    {
        EXCEPTION("A node index other than 0 or 1 has been requested for a Vessel Segment.");
    }

    if (mVessel.lock() != NULL)
    {
        mVessel.lock()->UpdateNodes();
    }
}

template<unsigned DIM>
template<typename T> void CaVesselSegment<DIM>::SetData(const std::string& variableName, T value)
{
    mDataContainer.SetData(variableName, value);
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetDataContainer(const VasculatureData& rDataContainer)
{
    mDataContainer = rDataContainer;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetLabel(const std::string& rLabel)
{
    mLabel = rLabel;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetRadius(double radius)
{
    mRadius = radius;
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVesselSegment<DIM>::Shared()
{
    boost::shared_ptr<CaVesselSegment<DIM> > pSegment = this->shared_from_this();
    return pSegment;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetFlowProperties(const SegmentFlowProperties& rFlowProperties)
{
    mpFlowProperties = boost::shared_ptr<SegmentFlowProperties>(new SegmentFlowProperties(rFlowProperties));
}

// Explicit instantiation
template class CaVesselSegment<2> ;
template class CaVesselSegment<3> ;

template bool CaVesselSegment<2>::GetData<bool>(const std::string& variableName) const;
template double CaVesselSegment<2>::GetData<double>(const std::string& variableName) const;
template unsigned CaVesselSegment<2>::GetData<unsigned>(const std::string& variableName) const;
template std::vector<double> CaVesselSegment<2>::GetData<std::vector<double> >(const std::string& variableName) const;
template void CaVesselSegment<2>::SetData(const std::string& variableName, bool value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, double value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, unsigned value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, std::vector<double> value);

template bool CaVesselSegment<3>::GetData<bool>(const std::string& variableName) const;
template double CaVesselSegment<3>::GetData<double>(const std::string& variableName) const;
template unsigned CaVesselSegment<3>::GetData<unsigned>(const std::string& variableName) const;
template std::vector<double> CaVesselSegment<3>::GetData<std::vector<double> >(const std::string& variableName) const;
template void CaVesselSegment<3>::SetData(const std::string& variableName, bool value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, double value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, unsigned value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, std::vector<double> value);
