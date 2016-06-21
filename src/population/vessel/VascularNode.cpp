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
#include "Exception.hpp"
#include "VascularNode.hpp"

template<unsigned DIM>
VascularNode<DIM>::VascularNode(double v1, double v2, double v3) :
        mLocation(ChastePoint<DIM>(v1 ,v2, v3)),
        mOutputData(),
        mId(0),
        mSegments(std::vector<boost::weak_ptr<VesselSegment<DIM> > >()),
        mRadius(10.0),
        mpFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false),
        mReferenceLength(1.0*unit::microns)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(c_vector<double, DIM> location) :
        mLocation(ChastePoint<DIM>(location)),
        mOutputData(),
        mId(0),
        mSegments(std::vector<boost::weak_ptr<VesselSegment<DIM> > >()),
        mRadius(10.0),
        mpFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false),
        mReferenceLength(1.0*unit::microns)
{
}

template<unsigned DIM>
VascularNode<DIM>::VascularNode(const VascularNode<DIM>& rExistingNode) :
        boost::enable_shared_from_this<VascularNode<DIM> >(),
        mLocation(rExistingNode.rGetLocation()),
        mOutputData(),
        mId(0),
        mSegments(std::vector<boost::weak_ptr<VesselSegment<DIM> > >()),
        mRadius(rExistingNode.GetRadius()),
        mpFlowProperties(boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties())),
        mIsMigrating(false),
        mReferenceLength(1.0*unit::microns)
{
    mOutputData = rExistingNode.GetOutputData();
    SetFlowProperties(*(rExistingNode.GetFlowProperties()));
    mIsMigrating = rExistingNode.IsMigrating();
}


template<unsigned DIM>
VascularNode<DIM>::~VascularNode()
{
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
void VascularNode<DIM>::AddSegment(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment)
{
    // Vessel segments can only be attached to a node once. Note use of lock to get shared_ptr from
    // weak_ptr.
    for (unsigned idx = 0; idx < mSegments.size(); idx++)
    {
        if (mSegments[idx].lock() == pVesselSegment)
        {
            EXCEPTION("This segment is already attached to this node.");
        }
    }
    mSegments.push_back(boost::weak_ptr<VesselSegment<DIM> >(pVesselSegment));
}

template<unsigned DIM>
double VascularNode<DIM>::GetDistance(const c_vector<double, DIM>& rLocation) const
{
    return norm_2(rLocation - mLocation.rGetLocation());
}

template<unsigned DIM>
units::quantity<unit::length> VascularNode<DIM>::GetDimensionalDistance(const c_vector<double, DIM>& rLocation) const
{
    return GetDistance(rLocation) * mReferenceLength;
}

template<unsigned DIM>
boost::shared_ptr<NodeFlowProperties> VascularNode<DIM>::GetFlowProperties() const
{
    return mpFlowProperties;
}

template<unsigned DIM>
unsigned VascularNode<DIM>::GetId() const
{
    return mId;
}

template<unsigned DIM>
const c_vector<double, DIM>& VascularNode<DIM>::rGetLocation() const
{
    return mLocation.rGetLocation();
}

template<unsigned DIM>
const c_vector<double, DIM>& VascularNode<DIM>::rGetDimensionalLocation() const
{
    return mLocation.rGetLocation()*(mReferenceLength/unit::metres);
}

template<unsigned DIM>
unsigned VascularNode<DIM>::GetNumberOfSegments() const
{
    return mSegments.size();
}

template<unsigned DIM>
double VascularNode<DIM>::GetOutputDataValue(const std::string& rKey) const
{
    std::map<std::string,double>::const_iterator it = mOutputData.find(rKey);
    if (it != mOutputData.end())
    {
        return it->second;
    }
    else
    {
        EXCEPTION("Requested output data key not found");
    }
}

template<unsigned DIM>
std::map<std::string, double> VascularNode<DIM>::GetOutputData() const
{
    return mOutputData;
}

template<unsigned DIM>
std::vector<std::string> VascularNode<DIM>::GetOutputDataKeys() const
{
    std::vector<std::string> keys;
    for(std::map<std::string,double>::const_iterator it = mOutputData.begin(); it != mOutputData.end(); ++it)
    {
        keys.push_back(it->first);
    }
    return keys;
}

template<unsigned DIM>
units::quantity<unit::length> VascularNode<DIM>::GetDimensionalRadius() const
{
    return mRadius * mReferenceLength;
}

template<unsigned DIM>
double VascularNode<DIM>::GetRadius() const
{
    return mRadius;
}

template<unsigned DIM>
units::quantity<unit::length> VascularNode<DIM>::GetReferenceLength() const
{
    return mReferenceLength;
}

template<unsigned DIM>
std::pair<double, std::string> VascularNode<DIM>::GetReferenceLengthValueAndUnit() const
{
    return LengthQuantityToValueUnitPair(mReferenceLength);
}

template<unsigned DIM>
boost::shared_ptr<VesselSegment<DIM> > VascularNode<DIM>::GetSegment(unsigned index) const
{
    if(index >= mSegments.size())
    {
        EXCEPTION("Requested segment index out of range");
    }
    else
    {
        // Convert to shared ptr from weak ptr
        return mSegments[index].lock();
    }
}

template<unsigned DIM>
std::vector<boost::shared_ptr<VesselSegment<DIM> > > VascularNode<DIM>::GetSegments() const
{
    // Need to do it this way because of weak pointers, can't just return mSegments
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments(mSegments.size());
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        segments[idx] = mSegments[idx].lock();
    }
    return segments;
}

template<unsigned DIM>
bool VascularNode<DIM>::IsAttachedTo(const boost::shared_ptr<VesselSegment<DIM> > pSegment) const
{
    // Need to get shared ptr from current node to allow for comparison
    return (pSegment->GetNode(0) == this->shared_from_this() || pSegment->GetNode(1) == this->shared_from_this());
}

template<unsigned DIM>
bool VascularNode<DIM>::IsCoincident(const c_vector<double, DIM>& rLocation) const
{
    bool returned_value = true;
    for (unsigned dim=0; dim<DIM; dim++)
    {
        if (rLocation[dim] != mLocation[dim])
        {
            returned_value = false;
            break;
        }
    }
    return returned_value;
}

template<unsigned DIM>
bool VascularNode<DIM>::IsMigrating() const
{
    return mIsMigrating;
}

template<unsigned DIM>
void VascularNode<DIM>::RemoveSegment(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment)
{
    // Need to do it this way due to weak pointer use
    for (unsigned idx = 0; idx < mSegments.size(); idx++)
    {
        if (mSegments[idx].lock() == pVesselSegment)
        {
            mSegments.erase(mSegments.begin() + idx);
            break;
        }
    }
}

template<unsigned DIM>
void VascularNode<DIM>::SetFlowProperties(const NodeFlowProperties& rFlowProperties)
{
    mpFlowProperties = boost::shared_ptr<NodeFlowProperties>(new NodeFlowProperties(rFlowProperties));
}

template<unsigned DIM>
void VascularNode<DIM>::SetLocation(const c_vector<double, DIM>& location)
{
    mLocation = ChastePoint<DIM>(location);
}

template<unsigned DIM>
void VascularNode<DIM>::SetLocation(double x, double y, double z)
{
    mLocation = ChastePoint<DIM>(x,y,z);
}

template<unsigned DIM>
void VascularNode<DIM>::SetOutputData(const std::string& rKey, double value)
{
    mOutputData[rKey] = value;
}

template<unsigned DIM>
void VascularNode<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void VascularNode<DIM>::SetIsMigrating(bool isMigrating)
{
    mIsMigrating = isMigrating;
}

template<unsigned DIM>
void VascularNode<DIM>::SetRadius(double radius)
{
    mRadius = radius;
}

template<unsigned DIM>
void VascularNode<DIM>::UpdateOutputData()
{
    std::map<std::string, double> flow_data = mpFlowProperties->GetVtkData();
    mOutputData.insert(flow_data.begin(), flow_data.end());
    mOutputData["Node Id"] = double(GetId());
    mOutputData["Dimensionless Node Radius"] = GetRadius();
    mOutputData["Dimensional Node Radius: " + LengthQuantityToString(mReferenceLength)] = GetDimensionalRadius();
    mOutputData["Node Is Migrating"] = double(IsMigrating());
}

// Explicit instantiation
template class VascularNode<2> ;
template class VascularNode<3> ;
