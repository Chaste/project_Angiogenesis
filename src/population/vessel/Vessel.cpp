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
#include "Exception.hpp"
#include "UblasIncludes.hpp"
#include "SimulationTime.hpp"

#include "Vessel.hpp"

template<unsigned DIM>
Vessel<DIM>::Vessel(boost::shared_ptr<VesselSegment<DIM> > pSegment) :
        mSegments(std::vector<boost::shared_ptr<VesselSegment<DIM> > >()),
        mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
        mNodesUpToDate(false),
        mDataContainer(),
        mId(0),
        mLabel(),
        mUndergoingRegression(false),
        mRemoveViaRegression(false),
        mRegressionTime(DBL_MAX)
{
    mSegments.push_back(pSegment);
}

template<unsigned DIM>
Vessel<DIM>::Vessel(std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments) :
        mSegments(segments),
        mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
        mNodesUpToDate(false),
        mDataContainer(),
        mId(0),
        mLabel(),
        mUndergoingRegression(false),
        mRemoveViaRegression(false),
        mRegressionTime(DBL_MAX)
{
    if (segments.size() > 1)
    {
        for (unsigned i = 1; i < mSegments.size(); i++)
        {
            if (!mSegments[i]->IsConnectedTo(mSegments[i - 1]))
            {
                EXCEPTION("Input vessel segments are not attached in the correct order.");
            }
        }

        for (unsigned i = 0; i < mSegments.size(); i++)
        {
            for (unsigned j = 0; j < mSegments.size(); j++)
            {
                if (i != j && i != j - 1 && i != j + 1)
                {
                    if (mSegments[i]->IsConnectedTo(mSegments[j]))
                    {
                        EXCEPTION("Input vessel segments are not correctly connected.");
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
Vessel<DIM>::Vessel(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes) :
        mSegments(std::vector<boost::shared_ptr<VesselSegment<DIM> > >()),
        mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
        mNodesUpToDate(false),
        mDataContainer(),
        mId(0),
        mLabel(),
        mUndergoingRegression(false),
        mRemoveViaRegression(false),
        mRegressionTime(DBL_MAX)
{

    if (nodes.size() < 2)
    {
        EXCEPTION("Insufficient number of nodes to define a segment.");
    }
    else
    {
        for (unsigned i = 1; i < nodes.size(); i++)
        {
            mSegments.push_back(VesselSegment<DIM>::Create(nodes[i-1], nodes[i]));
        }
    }
}

template<unsigned DIM>
Vessel<DIM>::Vessel(boost::shared_ptr<VascularNode<DIM> > pStartNode, boost::shared_ptr<VascularNode<DIM> > pEndNode)
    :        mSegments(std::vector<boost::shared_ptr<VesselSegment<DIM> > >()),
             mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
             mNodesUpToDate(false),
             mDataContainer(),
             mId(0),
             mLabel(),
             mUndergoingRegression(false),
             mRemoveViaRegression(false),
             mRegressionTime(DBL_MAX)
{
    mSegments.push_back(VesselSegment<DIM>::Create(pStartNode, pEndNode));
}

template<unsigned DIM>
Vessel<DIM>::~Vessel()
{
}

template<unsigned DIM>
boost::shared_ptr<Vessel<DIM> > Vessel<DIM>::Create(boost::shared_ptr<VesselSegment<DIM> > pSegment)
{
    boost::shared_ptr<Vessel<DIM> > pSelf(new Vessel<DIM>(pSegment));

    // Add the vessel to the segment
    pSegment->AddVessel(pSelf->shared_from_this());
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<Vessel<DIM> > Vessel<DIM>::Create(std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments)
{
    boost::shared_ptr<Vessel<DIM> > pSelf(new Vessel<DIM>(segments));

    // Add the vessel to the segments
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->AddVessel(pSelf->shared_from_this());
    }
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<Vessel<DIM> > Vessel<DIM>::Create(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes)
{
    boost::shared_ptr<Vessel<DIM> > pSelf(new Vessel<DIM>(nodes));

    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = pSelf->GetSegments();

    // Add the vessel to the new segments
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->AddVessel(pSelf->shared_from_this());
    }
    return pSelf;
}

template<unsigned DIM>
boost::shared_ptr<Vessel<DIM> > Vessel<DIM>::Create(boost::shared_ptr<VascularNode<DIM> > pStartNode, boost::shared_ptr<VascularNode<DIM> > pEndNode)
{
    boost::shared_ptr<Vessel<DIM> > pSelf(new Vessel<DIM>(pStartNode, pEndNode));

    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = pSelf->GetSegments();

    // Add the vessel to the new segment
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->AddVessel(pSelf->shared_from_this());
    }
    return pSelf;
}

template<unsigned DIM>
void Vessel<DIM>::AddSegment(boost::shared_ptr<VesselSegment<DIM> > pSegment)
{

    if(mSegments.size() == 1)
    {
        pSegment->AddVessel(Shared());
        if(pSegment->GetNode(0) == mSegments[0]->GetNode(0))
        {
            mSegments.insert(mSegments.begin(), pSegment);
        }
        else if(pSegment->GetNode(1) == mSegments[0]->GetNode(0))
        {
            mSegments.insert(mSegments.begin(), pSegment);
        }
        else if(pSegment->GetNode(0) == mSegments[0]->GetNode(1))
        {
            mSegments.push_back(pSegment);
        }
        else if(pSegment->GetNode(1) == mSegments[0]->GetNode(1))
        {
            mSegments.push_back(pSegment);
        }
        else
        {
            EXCEPTION("Input vessel segment does not coincide with any end of the vessel.");
        }
    }
    else
    {
        if (pSegment->IsConnectedTo(mSegments.back()))
        // Append to end of vessel
        {
            pSegment->AddVessel(Shared());
            mSegments.push_back(pSegment);
        }
        else if (pSegment->IsConnectedTo(mSegments.front()))
        // Insert at the start of the vessel
        {
            pSegment->AddVessel(Shared());
            mSegments.insert(mSegments.begin(), pSegment);
        }
        else
        {
            EXCEPTION("Input vessel segment does not coincide with any end of the multi-segment vessel.");
        }
    }

    mNodesUpToDate = false;
}

template<unsigned DIM>
void Vessel<DIM>::AddSegments(std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments)
{

    if (segments.front()->IsConnectedTo(mSegments.back()))
    {
        mSegments.insert(mSegments.end(), segments.begin(), segments.end());
    }
    else if (segments.back()->IsConnectedTo(mSegments.front()))
    {
        segments.insert(segments.end(), mSegments.begin(), mSegments.end());
        mSegments = segments;
    }
    else if (segments.front()->IsConnectedTo(mSegments.front()))
    {
        std::reverse(segments.begin(), segments.end());
        segments.insert(segments.end(), mSegments.begin(), mSegments.end());
        mSegments = segments;
    }
    else if (segments.back()->IsConnectedTo(mSegments.back()))
    {
        std::reverse(segments.begin(), segments.end());
        mSegments.insert(mSegments.end(), segments.begin(), segments.end());
    }
    else
    {
        EXCEPTION("Input vessel segments do not coincide with any end of the vessel.");
    }

    for (unsigned i = 1; i < mSegments.size(); i++)
    {
        if (!mSegments[i]->IsConnectedTo(mSegments[i - 1]))
        {
            EXCEPTION("Input vessel segments are not attached in the correct order.");
        }
    }

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        for (unsigned j = 0; j < mSegments.size(); j++)
        {
            if (i != j && i != j - 1 && i != j + 1)
            {
                if (mSegments[i]->IsConnectedTo(mSegments[j]))
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

    mNodesUpToDate = false;

}

template<unsigned DIM>
void Vessel<DIM>::CopyDataFromExistingVessel(boost::shared_ptr<Vessel<DIM> > pTargetVessel)
{
    mDataContainer.SetMap(pTargetVessel->rGetDataContainer().GetMap());
}

template<unsigned DIM>
template<typename T> T Vessel<DIM>::GetData(const std::string& variableName)
{
    return mDataContainer.GetData<T>(variableName);
}

template<unsigned DIM>
double Vessel<DIM>::GetClosestEndNodeDistance(c_vector<double, DIM> location)
{
    double distance_1 = this->GetStartNode()->GetDistance(location);
    double distance_2 = this->GetEndNode()->GetDistance(location);
    if(distance_1 > distance_2)
    {
        return distance_2;
    }
    else
    {
        return distance_1;
    }
}

template<unsigned DIM>
double Vessel<DIM>::GetDistance(c_vector<double, DIM> location) const
{
    // Get the distance to the nearest segment in the vessel
    double nearest_distance = DBL_MAX;
    for(unsigned idx=0; idx<mSegments.size(); idx++)
    {
        double seg_distance = mSegments[idx]->GetDistance(location);
        if(seg_distance < nearest_distance)
        {
            nearest_distance = seg_distance;
        }
    }
    return nearest_distance;
}

template<unsigned DIM>
const VasculatureData& Vessel<DIM>::rGetDataContainer() const
{
    return mDataContainer;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<Vessel<DIM> > > Vessel<DIM>::GetConnectedVessels()
{
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > start_segments = GetStartNode()->GetVesselSegments();
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > end_segments = GetStartNode()->GetVesselSegments();

    std::vector<boost::shared_ptr<Vessel<DIM> > > connected;

    for(unsigned idx=0; idx<start_segments.size();idx++)
    {
        if(start_segments[idx]->GetVessel() != Shared())
        {
            connected.push_back(start_segments[idx]->GetVessel());
        }
    }
    for(unsigned idx=0; idx<end_segments.size();idx++)
    {
        if(end_segments[idx]->GetVessel() != Shared())
        {
            connected.push_back(end_segments[idx]->GetVessel());
        }
    }
    return connected;
}

template<unsigned DIM>
std::vector<std::string> Vessel<DIM>::GetDataKeys(bool castable_to_double) const
{
    return mDataContainer.GetKeys(castable_to_double);
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > Vessel<DIM>::GetEndNode()
{
    if (!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes.back();
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > Vessel<DIM>::GetNodeAtOppositeEnd(
        boost::shared_ptr<VascularNode<DIM> > pQueryNode)
{
    if (!mNodesUpToDate)
    {
        UpdateNodes();
    }

    if (pQueryNode == GetStartNode())
    {
        return GetEndNode();
    }
    else if (pQueryNode == GetEndNode())
    {
        return GetStartNode();
    }
    else
    {
        EXCEPTION("Query node is not at either end of the vessel.");
    }
}

template<unsigned DIM>
unsigned Vessel<DIM>::GetId() const
{
    return mId;
}

template<unsigned DIM>
const std::string& Vessel<DIM>::rGetLabel() const
{
    return mLabel;
}

template<unsigned DIM>
double Vessel<DIM>::GetLength() const
{
    double length = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        length += mSegments[i]->GetLength();
    }

    return length;
}

template<unsigned DIM>
double Vessel<DIM>::GetRadius() const
{
    double radius = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        radius += mSegments[i]->GetRadius();
    }

    return radius / (double(mSegments.size()));
}

template<unsigned DIM>
double Vessel<DIM>::GetHaematocrit() const
{
    double haematocrit = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        haematocrit += mSegments[i]->GetFlowProperties()->GetHaematocrit();
    }

    return haematocrit / double(mSegments.size());
}

template<unsigned DIM>
double Vessel<DIM>::GetFlowRate() const
{
    double flow_rate = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        flow_rate += mSegments[i]->GetFlowProperties()->GetFlowRate();
    }

    return flow_rate / double(mSegments.size());
}

template<unsigned DIM>
double Vessel<DIM>::GetImpedance() const
{
    double impedance = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        impedance += mSegments[i]->GetFlowProperties()->GetImpedance();
    }

    return impedance;
}

template<unsigned DIM>
double Vessel<DIM>::GetViscosity() const
{
    double viscosity = 0.0;

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        viscosity += mSegments[i]->GetFlowProperties()->GetViscosity();
    }

    return viscosity / double(mSegments.size());
}

template<unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > Vessel<DIM>::GetNodes()
{
    if (!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes;
}

template<unsigned DIM>
unsigned Vessel<DIM>::GetNumberOfNodes()
{
    if (!mNodesUpToDate)
    {
        UpdateNodes();
    }
    return GetNumberOfSegments() + 1;
}

template<unsigned DIM>
unsigned Vessel<DIM>::GetNumberOfSegments()
{
    return mSegments.size();
}

template<unsigned DIM>
boost::shared_ptr<VesselSegment<DIM> > Vessel<DIM>::GetSegment(unsigned i)
{
    if (i < mSegments.size())
    {
        return mSegments[i];
    }
    else
    {
        EXCEPTION("Requested segment index exceeds number of segments.");
    }
}

template<unsigned DIM>
std::vector<boost::shared_ptr<VesselSegment<DIM> > > Vessel<DIM>::GetSegments()
{
    return mSegments;
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > Vessel<DIM>::GetStartNode()
{
    if (!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes.front();
}

template<unsigned DIM>
std::map<std::string, double> Vessel<DIM>::GetVtkData() const
{
    std::map<std::string, double> vtk_data;
    vtk_data["Vessel Id"] = double(GetId());
    vtk_data["Vessel Radius"] = GetRadius();
    vtk_data["Vessel Impedance"] = GetImpedance();
    vtk_data["Vessel Haematocrit"] = GetHaematocrit();
    vtk_data["Vessel Flow Rate"] = GetFlowRate();
    vtk_data["Vessel Absolute Flow Rate"] = std::abs(GetFlowRate());
    vtk_data["Vessel Viscosity"] = GetViscosity();
    return vtk_data;
}

template<unsigned DIM>
bool Vessel<DIM>::HasDataKey(const std::string& rKey) const
{
    return mDataContainer.HasKey(rKey);
}

template<unsigned DIM>
bool Vessel<DIM>::IsConnectedTo(boost::shared_ptr<Vessel<DIM> > pOtherVessel)
{
    if (GetStartNode() == pOtherVessel->GetStartNode() || GetEndNode() == pOtherVessel->GetStartNode()
            || GetStartNode() == pOtherVessel->GetEndNode() || GetEndNode() == pOtherVessel->GetEndNode())
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > Vessel<DIM>::DivideSegment(ChastePoint<DIM> location)
{
    // Identify segment
    boost::shared_ptr<VesselSegment<DIM> > pVesselSegment;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        if (mSegments[i]->GetDistance(location) <= 1e-6)
        {
            pVesselSegment = mSegments[i];
            if (pVesselSegment->GetNode(0)->IsCoincident(location))
            {
                return pVesselSegment->GetNode(0);
            }
            if (pVesselSegment->GetNode(1)->IsCoincident(location))
            {
                return pVesselSegment->GetNode(1);
            }

        }
    }
    if (!pVesselSegment)
    {
        EXCEPTION("Specified location is not on a segment in this vessel.");
    }
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        for (unsigned j = 0; j < mSegments.size(); j++)
        {
            if (i != j && i != j - 1 && i != j + 1)
            {
                if (mSegments[i]->IsConnectedTo(mSegments[j]))
                {
                    EXCEPTION("Input vessel segments are not correctly connected.");
                }
            }
        }
    }

    // The node's data is averaged from the original segments's nodes
    // Get the closest node
    double distance0 = pVesselSegment->GetNode(0)->GetDistance(location);
    double distance1 = pVesselSegment->GetNode(1)->GetDistance(location);
    unsigned closest_index;

    if (distance0 <= distance1)
    {
        closest_index = 0;
    }
    else
    {
        closest_index = 1;
    }

    // Make a copy of the closest node
    boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(*pVesselSegment->GetNode(closest_index));
    p_new_node->SetLocation(location);

    // Make two new segments
    boost::shared_ptr<VesselSegment<DIM> > p_new_segment0 =
          VesselSegment<DIM>::Create(pVesselSegment->GetNode(0), p_new_node);

    boost::shared_ptr<VesselSegment<DIM> > p_new_segment1 =
            VesselSegment<DIM>::Create(p_new_node, pVesselSegment->GetNode(1));

    p_new_segment0->CopyDataFromExistingSegment(pVesselSegment);
    p_new_segment1->CopyDataFromExistingSegment(pVesselSegment);

    // Add new segments, ensuring they are correctly ordered
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > newSegments;
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it = std::find(mSegments.begin(),
                                                                                             mSegments.end(),
                                                                                             pVesselSegment);

    if (it == mSegments.end())
    {
        EXCEPTION("Vessel segment is not contained inside vessel.");
    }

    if (mSegments.size() == 1)
    {
        newSegments.push_back(p_new_segment0);
        newSegments.push_back(p_new_segment1);
    }
    else if(it == mSegments.begin())
    {
        if((*(it + 1))->IsConnectedTo(p_new_segment1))
        {
            newSegments.push_back(p_new_segment0);
            newSegments.push_back(p_new_segment1);
        }
        else
        {
            newSegments.push_back(p_new_segment1);
            newSegments.push_back(p_new_segment0);
        }
    }
    else if ((*(it - 1))->IsConnectedTo(p_new_segment0))
    {
        newSegments.push_back(p_new_segment0);
        newSegments.push_back(p_new_segment1);
    }
    else
    {
        newSegments.push_back(p_new_segment1);
        newSegments.push_back(p_new_segment0);
    }
    mSegments.insert(it, newSegments.begin(), newSegments.end());

    // Remove old segment
    it = std::find(mSegments.begin(), mSegments.end(), pVesselSegment);

    if (it != mSegments.end())
    {
        mSegments.erase(it);
        pVesselSegment->Remove();
    }
    else
    {
        EXCEPTION("Vessel segment is not contained inside vessel.");
    }
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        for (unsigned j = 0; j < mSegments.size(); j++)
        {
            if (i != j && i != j - 1 && i != j + 1)
            {
                if (mSegments[i]->IsConnectedTo(mSegments[j]))
                {
                    EXCEPTION("Input vessel segments are not correctly connected.");
                }
            }
        }
    }

    // Add the vessel to the segments
    for (unsigned i = 0; i < newSegments.size(); i++)
    {
        newSegments[i]->AddVessel(Shared());
    }

    mNodesUpToDate = false;
    return p_new_node;
}

template<unsigned DIM>
void Vessel<DIM>::Remove()
{
    // Detach all segments from their nodes
    for (unsigned idx = 0; idx < mSegments.size(); idx++)
    {
        mSegments[idx]->Remove();
    }
    mNodesUpToDate = false;
    mSegments = std::vector<boost::shared_ptr<VesselSegment<DIM> > >();
}

template<unsigned DIM>
void Vessel<DIM>::RemoveSegments(SegmentLocation::Value location)
{

    if (mSegments.size() == 1)
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

    mNodesUpToDate = false;
}

template<unsigned DIM>
template<typename T> void Vessel<DIM>::SetData(const std::string& variableName, T value)
{
    mDataContainer.SetData(variableName, value);
}

template<unsigned DIM>
void Vessel<DIM>::SetDataContainer(const VasculatureData& rDataContainer)
{
    mDataContainer = rDataContainer;
}

template<unsigned DIM>
void Vessel<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void Vessel<DIM>::SetLabel(const std::string& label)
{
    mLabel = label;
}

template<unsigned DIM>
void Vessel<DIM>::SetRadius(double radius)
{
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->SetRadius(radius);
    }
}

template<unsigned DIM>
void Vessel<DIM>::SetHaematocrit(double haematocrit)
{
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetHaematocrit(haematocrit);
    }
}

template<unsigned DIM>
void Vessel<DIM>::SetFlowRate(double flowRate)
{
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetFlowRate(flowRate);
    }
}

template<unsigned DIM>
boost::shared_ptr<Vessel<DIM> > Vessel<DIM>::Shared()
{
    return this->shared_from_this();
}

template<unsigned DIM>
void Vessel<DIM>::UpdateNodes()
{
    mNodes = std::vector<boost::shared_ptr<VascularNode<DIM> > >();

    if (mSegments.size() == 1)
    {
        mNodes.push_back(mSegments[0]->GetNode(0));
        mNodes.push_back(mSegments[0]->GetNode(1));
    }

    // Add the start and end nodes of the first segment and then
    // the end nodes of every other segment
    else
    {
        if (mSegments[1]->HasNode(mSegments[0]->GetNode(1)))
        {
            mNodes.push_back(mSegments[0]->GetNode(0));
            mNodes.push_back(mSegments[0]->GetNode(1));
        }
        else if (mSegments[1]->HasNode(mSegments[0]->GetNode(1)))
        {
            mNodes.push_back(mSegments[0]->GetNode(1));
            mNodes.push_back(mSegments[0]->GetNode(0));

        }

        for (unsigned idx = 1; idx < mSegments.size(); idx++)
        {
            if (mNodes[idx] == mSegments[idx]->GetNode(0))
            {
                mNodes.push_back(mSegments[idx]->GetNode(1));
            }
            else
            {
                mNodes.push_back(mSegments[idx]->GetNode(0));
            }
        }
    }
    mNodesUpToDate = true;
}


template<unsigned DIM>
void Vessel<DIM>::SetTimeUntilRegression(double time)
{
    assert(!mRemoveViaRegression);
    assert(time >= 0);

    if (HasRegressionTimerStarted())
    {
        EXCEPTION("SetTimeUntilRegression(time) called when already undergoing regression");
    }

    mUndergoingRegression = true;
    mRegressionTime = SimulationTime::Instance()->GetTime() + time;
}

template<unsigned DIM>
bool Vessel<DIM>::HasRegressionTimerStarted()
{
    return mUndergoingRegression;
}

template<unsigned DIM>
void Vessel<DIM>::ResetRegressionTimer()
{
    mRegressionTime = DBL_MAX;
    mUndergoingRegression = false;
    mRemoveViaRegression = false;
}

template<unsigned DIM>
bool Vessel<DIM>::VesselHasRegressed()
{
    if (SimulationTime::Instance()->GetTime() >= mRegressionTime)
    {
        assert(mUndergoingRegression);
        mRemoveViaRegression = true;
        return mRemoveViaRegression;
    }
    else
    {
        mRemoveViaRegression = false;
        return mRemoveViaRegression;
    }
}



// Explicit instantiation
template class Vessel<2> ;
template class Vessel<3> ;

template bool Vessel<2>::GetData<bool>(const std::string& variableName);
template double Vessel<2>::GetData<double>(const std::string& variableName);
template unsigned Vessel<2>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> Vessel<2>::GetData<std::vector<double> >(const std::string& variableName);
template void Vessel<2>::SetData(const std::string& variableName, bool value);
template void Vessel<2>::SetData(const std::string& variableName, double value);
template void Vessel<2>::SetData(const std::string& variableName, unsigned value);
template void Vessel<2>::SetData(const std::string& variableName, std::vector<double> value);

template bool Vessel<3>::GetData<bool>(const std::string& variableName);
template double Vessel<3>::GetData<double>(const std::string& variableName);
template unsigned Vessel<3>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> Vessel<3>::GetData<std::vector<double> >(const std::string& variableName);
template void Vessel<3>::SetData(const std::string& variableName, bool value);
template void Vessel<3>::SetData(const std::string& variableName, double value);
template void Vessel<3>::SetData(const std::string& variableName, unsigned value);
template void Vessel<3>::SetData(const std::string& variableName, std::vector<double> value);
