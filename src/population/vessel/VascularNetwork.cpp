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

#include <iostream>
#include <math.h>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the boost deprecated warning for now (gcc4.3)
#include <boost/graph/graphviz.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#ifdef CHASTE_VTK
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#endif // CHASTE_VTK
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "SegmentFlowProperties.hpp"

#include "VascularNetwork.hpp"

template <unsigned DIM>
VascularNetwork<DIM>::VascularNetwork()
: mVessels(std::vector<boost::shared_ptr<Vessel<DIM> > >()),
  mSegments(std::vector<boost::shared_ptr<VesselSegment<DIM> > >()),
  mSegmentsUpToDate(false),
  mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
  mNodesUpToDate(false),
  mVesselNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
  mVesselNodesUpToDate(false),
  mDataContainer()
{

}

template <unsigned DIM>
VascularNetwork<DIM>::~VascularNetwork()
{

}

template <unsigned DIM>
boost::shared_ptr<VascularNetwork<DIM> > VascularNetwork<DIM>::Create()
{
    MAKE_PTR(VascularNetwork<DIM>, pSelf);
    return pSelf;
}

template <unsigned DIM>
void VascularNetwork<DIM>::AddVessel(boost::shared_ptr<Vessel<DIM> > pVessel)
{
    mVessels.push_back(pVessel);
    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void VascularNetwork<DIM>::AddVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels)
{
    mVessels.insert(mVessels.end(), vessels.begin(), vessels.end());
    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void VascularNetwork<DIM>::CopySegmentFlowProperties(unsigned index)
{
    boost::shared_ptr<SegmentFlowProperties> properties = GetVesselSegments()[index]->GetFlowProperties();
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetFlowProperties(*properties);
    }
}

template <unsigned DIM>
std::vector<boost::shared_ptr<Vessel<DIM> > > VascularNetwork<DIM>::CopyVessels()
{
    return CopyVessels(mVessels);
}

template <unsigned DIM>
std::vector<boost::shared_ptr<Vessel<DIM> > > VascularNetwork<DIM>::CopyVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels)
{
    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator vessel_iter;
    std::vector<boost::shared_ptr<Vessel<DIM> > > new_vessels;
    for(vessel_iter = vessels.begin(); vessel_iter != vessels.end(); vessel_iter++)
    {
        typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator segment_iter;
        std::vector<boost::shared_ptr<VesselSegment<DIM> > > new_segments;
        std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = (*vessel_iter)->GetSegments();
        for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
        {
            ChastePoint<DIM> node0_location = (*segment_iter)->GetNode(0)->GetLocation();
            ChastePoint<DIM> node1_location = (*segment_iter)->GetNode(1)->GetLocation();
            new_segments.push_back(VesselSegment<DIM>::Create(VascularNode<DIM>::Create(node0_location),
                    VascularNode<DIM>::Create(node1_location)));
        }
        new_vessels.push_back(Vessel<DIM>::Create(new_segments));
    }

    MergeCoincidentNodes(new_vessels);
    AddVessels(new_vessels);
    return new_vessels;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNetwork<DIM>::DivideVessel(boost::shared_ptr<Vessel<DIM> > pVessel, ChastePoint<DIM> location)
{
    boost::shared_ptr<VesselSegment<DIM> > p_segment;

    // If the divide location coincides with one of the end nodes don't divide and return that node
    if (pVessel->GetStartNode()->IsCoincident(location) || pVessel->GetEndNode()->IsCoincident(location))
    {

        if (pVessel->GetStartNode()->IsCoincident(location))
        {
            return pVessel->GetStartNode();
        }
        else
        {
            return pVessel->GetEndNode();
        }
    }
    else
    {
        bool locatedInsideVessel = false;
        std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = pVessel->GetSegments();
        for (unsigned idx = 0; idx < segments.size(); idx++)
        {
            if (segments[idx]->GetDistance(location) <= 1e-6)
            {
                locatedInsideVessel = true;
                p_segment = segments[idx];
                break;
            }
        }
        if(!locatedInsideVessel)
        {
            EXCEPTION("There is no segment at the requested division location.");
        }
    }
    boost::shared_ptr<VascularNode<DIM> > p_new_node = pVessel->DivideSegment(location); // network segments and nodes out of date

    // create two new vessels and assign them the old vessel's properties
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > start_segments;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > end_segments;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = pVessel->GetSegments();
    unsigned segment_index = segments.size()+1;
    for (unsigned idx = 0; idx < segments.size(); idx++)
    {
        start_segments.push_back(segments[idx]);
        if (segments[idx]->GetNode(1)->IsCoincident(location))
        {
            segment_index = idx;
            break;
        }
    }

    if (segment_index == segments.size()-1)
    {
        EXCEPTION("Vessel segment not found.");
    }
    for (unsigned idx = segment_index+1; idx < segments.size(); idx++)
    {
        end_segments.push_back(segments[idx]);
    }
    boost::shared_ptr<Vessel<DIM> > p_new_vessel1 = Vessel<DIM>::Create(start_segments);
    boost::shared_ptr<Vessel<DIM> > p_new_vessel2 = Vessel<DIM>::Create(end_segments);
    p_new_vessel1->CopyDataFromExistingVessel(pVessel);
    p_new_vessel2->CopyDataFromExistingVessel(pVessel);

    AddVessel(p_new_vessel1);
    AddVessel(p_new_vessel2);
    RemoveVessel(pVessel, false);

    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
    return p_new_node;
}

template <unsigned DIM>
void VascularNetwork<DIM>::ExtendVessel(boost::shared_ptr<Vessel<DIM> > pVessel,
                  boost::shared_ptr<VascularNode<DIM> > pEndNode,
                  boost::shared_ptr<VascularNode<DIM> > pNewNode)
{
    if(pVessel->GetStartNode() == pEndNode)
    {
        boost::shared_ptr<VesselSegment<DIM> > p_segment = VesselSegment<DIM>::Create(pNewNode, pEndNode);
        p_segment->SetFlowProperties(*(pEndNode->GetVesselSegment(0)->GetFlowProperties()));
        p_segment->SetRadius(pEndNode->GetVesselSegment(0)->GetRadius());
        pVessel->AddSegment(p_segment);
    }
    else
    {
        boost::shared_ptr<VesselSegment<DIM> > p_segment = VesselSegment<DIM>::Create(pEndNode, pNewNode);
        p_segment->SetFlowProperties(*(pEndNode->GetVesselSegment(0)->GetFlowProperties()));
        p_segment->SetRadius(pEndNode->GetVesselSegment(0)->GetRadius());
        pVessel->AddSegment(p_segment);
    }

    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
boost::shared_ptr<Vessel<DIM> > VascularNetwork<DIM>::FormSprout(ChastePoint<DIM> sproutBaseLocation, ChastePoint<DIM> sproutTipLocation)
{
    // locate vessel at which the location of the sprout base exists
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> nearest_segment = GetNearestSegment(sproutBaseLocation);
    if (nearest_segment.second > 1e-6)
    {
        EXCEPTION("No vessel located at sprout base.");
    }

    // divide vessel at location of sprout base
    boost::shared_ptr<VascularNode<DIM> > p_new_node = DivideVessel(nearest_segment.first->GetVessel(), sproutBaseLocation);

    // create new vessel
    boost::shared_ptr<VascularNode<DIM> > p_new_node_at_tip = VascularNode<DIM>::Create(p_new_node);
    p_new_node_at_tip->SetLocation(sproutTipLocation);
    p_new_node_at_tip->SetIsMigrating(true);
    boost::shared_ptr<VesselSegment<DIM> > p_new_segment = VesselSegment<DIM>::Create(p_new_node, p_new_node_at_tip);
    p_new_segment->CopyDataFromExistingSegment(nearest_segment.first);
    boost::shared_ptr<Vessel<DIM> > p_new_vessel = Vessel<DIM>::Create(p_new_segment);
    AddVessel(p_new_vessel);
    return p_new_vessel;
}

template <unsigned DIM>
std::vector<std::pair<double, double> > VascularNetwork<DIM>::GetExtents()
{
    double x_max = -DBL_MAX;
    double y_max = -DBL_MAX;
    double z_max = -DBL_MAX;

    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        ChastePoint<DIM> location = (*it)->GetLocation();
        if(location[0] > x_max)
        {
            x_max = location[0];
        }
        if(location[1] > y_max)
        {
            y_max = location[1];
        }
        if(DIM > 2)
        {
            if(location[2] > z_max)
            {
                z_max = location[2];
            }
        }
    }

    double x_min = x_max;
    double y_min = y_max;
    double z_min = z_max;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        ChastePoint<DIM> location = (*it)->GetLocation();
        if(location[0] < x_min)
        {
            x_min = location[0];
        }
        if(location[1] < y_min)
        {
            y_min = location[1];
        }
        if(DIM > 2)
        {
            if(location[2] < z_min)
            {
                z_min = location[2];
            }
        }
    }

    std::vector<std::pair<double, double> > container;
    container.push_back(std::pair<double, double>(x_min, x_max));
    container.push_back(std::pair<double, double>(y_min, y_max));
    container.push_back(std::pair<double, double>(z_min, z_max));
    return container;
}

template <unsigned DIM>
double VascularNetwork<DIM>::GetDistanceToNearestNode(const ChastePoint<DIM>& rLocation)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    boost::shared_ptr<VascularNode<DIM> > nearest_node;
    double min_distance = 1.e12;

    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {

        double node_distance = (*node_iter)->GetDistance(rLocation);
        if (node_distance < min_distance)
        {
            min_distance = node_distance;
            nearest_node = (*node_iter) ;
        }
    }

    return min_distance;
}

template <unsigned DIM>
std::vector<double> VascularNetwork<DIM>::GetInterCapillaryDistances()
{
    std::vector<double> distances;
    for(unsigned idx=0; idx<mVessels.size(); idx++)
    {
        double min_distance = 1.e6;
        for(unsigned jdx=0; jdx<mVessels.size(); jdx++)
        {
            if(mVessels[idx]!=mVessels[jdx])
            {
                double distance = mVessels[idx]->GetStartNode()->GetDistance(mVessels[jdx]->GetStartNode());
                if(distance < min_distance)
                {
                    min_distance = distance;
                }
            }
        }
        distances.push_back(min_distance);
    }
    return distances;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNetwork<DIM>::GetNearestNode(const ChastePoint<DIM>& rLocation)
{
    return GetNearestNode(rLocation.rGetLocation());
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNetwork<DIM>::GetNearestNode(boost::shared_ptr<VascularNode<DIM> > pInputNode)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    boost::shared_ptr<VascularNode<DIM> > nearest_node;
    double min_distance = 1.e12;

    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {
        if((*node_iter) != pInputNode)
        {
            double node_distance = (*node_iter)->GetDistance(pInputNode->GetLocationVector());
            if (node_distance < min_distance)
            {
                min_distance = node_distance;
                nearest_node = (*node_iter) ;
            }
        }
    }
    return nearest_node;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNetwork<DIM>::GetNearestNode(c_vector<double, DIM> location)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    boost::shared_ptr<VascularNode<DIM> > nearest_node;
    double min_distance = 1.e12;

    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {

        double node_distance = (*node_iter)->GetDistance(location);
        if (node_distance < min_distance)
        {
            min_distance = node_distance;
            nearest_node = (*node_iter) ;
        }
    }

    return nearest_node;
}

template <unsigned DIM>
std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> VascularNetwork<DIM>::GetNearestSegment(boost::shared_ptr<VesselSegment<DIM> > pSegment)
{
    boost::shared_ptr<VesselSegment<DIM> > nearest_segment;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    double min_distance = 1.e12;
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator segment_iter;
    for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
    {
        if(!pSegment->IsConnectedTo((*segment_iter)))
        {
            // Get the segment to segment distance (http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment())
            c_vector<double, DIM> u = (*segment_iter)->GetNode(1)->GetLocationVector() - (*segment_iter)->GetNode(0)->GetLocationVector();
            c_vector<double, DIM> v = pSegment->GetNode(1)->GetLocationVector() - pSegment->GetNode(0)->GetLocationVector();
            c_vector<double, DIM> w = (*segment_iter)->GetNode(0)->GetLocationVector() - pSegment->GetNode(0)->GetLocationVector();

            double a = inner_prod(u,u);
            double b = inner_prod(u,v);
            double c = inner_prod(v,v);
            double d = inner_prod(u,w);
            double e = inner_prod(v,w);

            double dv = a * c - b * b;
            double sc, sn, sd = dv;
            double tc, tn ,td = dv;

            if(dv < 1.e-12) // almost parallel segments
            {
                sn = 0.0;
                sd = 1.0;
                tn = e;
                td = c;
            }
            else // get the closest point on the equivalent infinite lines
            {
                sn = (b*e - c*d);
                tn = (a*e - b*d);
                if ( sn < 0.0)
                {
                    sn = 0.0;
                    tn = e;
                    td = c;
                }
                else if(sn > sd)
                {
                    sn =sd;
                    tn = e+ b;
                    td = c;
                }
            }

            if(tn < 0.0)
            {
                tn = 0.0;
                if(-d < 0.0)
                {
                    sn = 0.0;
                }
                else if(-d > a)
                {
                    sn = sd;
                }
                else
                {
                    sn = -d;
                    sd = a;
                }
            }
            else if(tn > td)
            {
                tn = td;
                if((-d + b) < 0.0)
                {
                    sn = 0.0;
                }
                else if((-d + b) > a)
                {
                    sn = sd;
                }
                else
                {
                    sn = (-d + b);
                    sd = a;
                }
            }

            sc = (std::abs(sn) < 1.e-12 ? 0.0 : sn/sd);
            tc = (std::abs(tn) < 1.e-12 ? 0.0 : tn/td);
            c_vector<double, DIM> dp = w + (sc * u) - (tc * v);

            double segment_distance = norm_2(dp);
            if (segment_distance < min_distance)
            {
                min_distance = segment_distance;
                nearest_segment = (*segment_iter) ;
            }
        }

    }
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> return_pair =
            std::pair<boost::shared_ptr<VesselSegment<DIM> >, double>(nearest_segment, min_distance);
    return return_pair;
}

template <unsigned DIM>
std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> VascularNetwork<DIM>::GetNearestSegment(boost::shared_ptr<VascularNode<DIM> > pNode, bool sameVessel)
{
    boost::shared_ptr<VesselSegment<DIM> > nearest_segment;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    double min_distance = 1.e12;
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator segment_iter;
    for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
    {
        double segment_distance = (*segment_iter)->GetDistance(pNode->GetLocationVector());

        if (segment_distance < min_distance && (*segment_iter)->GetNode(0) != pNode && (*segment_iter)->GetNode(1) != pNode)
        {
            if(sameVessel)
            {
                min_distance = segment_distance;
                nearest_segment = (*segment_iter) ;
            }
            else
            {
                bool same_vessel = false;
                std::vector<boost::shared_ptr<VesselSegment<DIM> > > node_segs = pNode->GetVesselSegments();
                for(unsigned idx=0;idx<node_segs.size();idx++)
                {
                    if(node_segs[idx]->GetVessel() == (*segment_iter)->GetVessel())
                    {
                        same_vessel = true;
                    }
                }
                if(!same_vessel)
                {
                    min_distance = segment_distance;
                    nearest_segment = (*segment_iter);
                }
            }
        }
    }
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> return_pair =
            std::pair<boost::shared_ptr<VesselSegment<DIM> >, double>(nearest_segment, min_distance);
    return return_pair;
}

template <unsigned DIM>
std::pair<boost::shared_ptr<VesselSegment<DIM> >, double>  VascularNetwork<DIM>::GetNearestSegment(const ChastePoint<DIM>& rLocation)
{
    return GetNearestSegment(rLocation.rGetLocation());
}

template <unsigned DIM>
std::pair<boost::shared_ptr<VesselSegment<DIM> >, double>  VascularNetwork<DIM>::GetNearestSegment(c_vector<double, DIM> location)
{
    boost::shared_ptr<VesselSegment<DIM> > nearest_segment;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    double min_distance = 1.e12;
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator segment_iter;
    for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
    {
        double segment_distance = (*segment_iter)->GetDistance(location);
        if (segment_distance < min_distance)
        {
            min_distance = segment_distance;
            nearest_segment = (*segment_iter) ;
        }
    }
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> return_pair =
            std::pair<boost::shared_ptr<VesselSegment<DIM> >, double>(nearest_segment, min_distance);
    return return_pair;
}

template <unsigned DIM>
boost::shared_ptr<Vessel<DIM> > VascularNetwork<DIM>::GetNearestVessel(const ChastePoint<DIM>& rLocation)
{
    return GetNearestSegment(rLocation).first->GetVessel();
}

template <unsigned DIM>
boost::shared_ptr<Vessel<DIM> > VascularNetwork<DIM>::GetNearestVessel(c_vector<double, DIM> location)
{
    return GetNearestSegment(location).first->GetVessel();
}

template <unsigned DIM>
double VascularNetwork<DIM>::GetTotalLength()
{
    double length = 0.0;
    for(unsigned idx=0; idx<mVessels.size(); idx++)
    {
        length += mVessels[idx]->GetLength();
    }
    return length;
}

template <unsigned DIM>
double VascularNetwork<DIM>::GetTotalVolume()
{
    double volume = 0.0;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();
    for(unsigned idx=0; idx< segments.size(); idx++)
    {
        volume += segments[idx]->GetLength() * segments[idx]->GetRadius() * segments[idx]->GetRadius() * M_PI;
    }
    return volume;

}

template <unsigned DIM>
double VascularNetwork<DIM>::GetTotalSurfaceArea()
{
    double area = 0.0;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();
    for(unsigned idx=0; idx< segments.size(); idx++)
    {
        area += segments[idx]->GetLength() * 2.0 * segments[idx]->GetRadius() * M_PI;
    }
    return area;
}

template <unsigned DIM>
double VascularNetwork<DIM>::GetAverageInterSegmentDistance()
{
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    // store segment midpoints
    std::vector<c_vector<double, DIM> > midpoints(segments.size());
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        midpoints[idx] = segments[idx]->GetMidPoint();
    }

    // get intersegment distances
    double av_dist = 0.0;
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        double min_dist = 1.e6;
        for(unsigned jdx=0; jdx<segments.size(); jdx++)
        {
            if(segments[idx] != segments[jdx] && segments[idx]->GetVessel() != segments[jdx]->GetVessel())
            {
                double dist = norm_2(midpoints[idx] - midpoints[jdx]);
                if(dist < min_dist)
                {
                    min_dist = dist;
                }
            }
        }
        av_dist += min_dist;
    }
    return av_dist / double(segments.size());
}

template <unsigned DIM>
double VascularNetwork<DIM>::GetAverageVesselLength()
{
    return GetTotalLength() / double(mVessels.size());
}

template <unsigned DIM>
std::vector<unsigned> VascularNetwork<DIM>::GetVesselLengthDistribution(double binSpacing, unsigned numberOfBins)
{
    std::vector<unsigned> bins(numberOfBins, 0);

    // populate the bins
    for(unsigned idx=0; idx<mVessels.size(); idx++)
    {
        unsigned bin_label = std::floor(mVessels[idx]->GetLength() / binSpacing);
        if(bin_label > numberOfBins)
        {
            bin_label = numberOfBins;
        }
        bins[bin_label]++;
    }
    return bins;
}

template <unsigned DIM>
void VascularNetwork<DIM>::RemoveShortVessels(double cutoff, bool endsOnly)
{
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels_to_remove;

    for(unsigned idx=0; idx<mVessels.size(); idx++)
    {
        if(mVessels[idx]->GetLength() < cutoff)
        {
            if(endsOnly && (mVessels[idx]->GetStartNode()->GetNumberOfSegments() == 1 || mVessels[idx]->GetEndNode()->GetNumberOfSegments() == 1 ))
            {
                vessels_to_remove.push_back(mVessels[idx]);
            }
            else if(!endsOnly)
            {
                vessels_to_remove.push_back(mVessels[idx]);
            }
        }
    }

    for(unsigned idx=0; idx<vessels_to_remove.size(); idx++)
    {
        RemoveVessel(vessels_to_remove[idx]);
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::MergeShortVessels(double cutoff)
{
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels_to_merge;
    for(unsigned idx=0; idx<mVessels.size(); idx++)
    {
        if(mVessels[idx]->GetLength() < cutoff)
        {
            vessels_to_merge.push_back(mVessels[idx]);
        }
    }

    // Get the nodes, remove the vessel, move the nodes together
    for(unsigned idx=0; idx<vessels_to_merge.size(); idx++)
    {
        vessels_to_merge[idx]->GetEndNode()->SetLocation(vessels_to_merge[idx]->GetStartNode()->GetLocationVector());
        RemoveVessel(vessels_to_merge[idx], true);
    }

    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
    MergeCoincidentNodes();
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::NumberOfNodesNearLocation(const ChastePoint<DIM>& rLocation, double radius)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    unsigned num_nodes = 0;

    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        if(nodes[idx]->GetDistance(rLocation) <= radius + 1.e-6)
        {
            num_nodes++;
        }
    }
    return num_nodes;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > VascularNetwork<DIM>::GetNodes()
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > VascularNetwork<DIM>::GetNode(unsigned index)
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes[index];
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetNumberOfNodes()
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes.size();
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetNumberOfVesselNodes()
{
    if(!mVesselNodesUpToDate)
    {
        UpdateVesselNodes();
    }

    return mVesselNodes.size();
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetNodeIndex(boost::shared_ptr<VascularNode<DIM> > node)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > vec_nodes = GetNodes();

    unsigned index = 0;

    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it;
    for(it = vec_nodes.begin(); it != vec_nodes.end(); it++)
    {
        if(*it == node)
        {
            return index;
        }
        index ++;
    }
    EXCEPTION("Node is not in the network.");

}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetNumberOfVessels()
{
    return mVessels.size();
}

template <unsigned DIM>
boost::shared_ptr<Vessel<DIM> > VascularNetwork<DIM>::GetVessel(unsigned index)
{
    return mVessels[index];
}

template <unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > VascularNetwork<DIM>::GetVesselEndNodes()
{

    if(!mVesselNodesUpToDate)
    {
        UpdateVesselNodes();
    }

    return mVesselNodes;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<Vessel<DIM> > > VascularNetwork<DIM>::GetVessels()
{
    return mVessels;
}

template <unsigned DIM>
std::vector<std::vector<unsigned> > VascularNetwork<DIM>::GetNodeNodeConnectivity()
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels = GetVessels();
    std::vector<std::vector<unsigned> > node_vessel_connectivity = GetNodeVesselConnectivity();

    std::vector<std::vector<unsigned> > connectivity;
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        std::vector<unsigned> node_indexes;
        boost::shared_ptr<VascularNode<DIM> > p_node = nodes[node_index];
        unsigned num_branches = node_vessel_connectivity[node_index].size();
        for (unsigned vessel_index = 0; vessel_index < num_branches; vessel_index++)
        {
            boost::shared_ptr<Vessel<DIM> > p_vessel = vessels[node_vessel_connectivity[node_index][vessel_index]];

            // Get the node at the other end of the vessel
            boost::shared_ptr<VascularNode<DIM> > p_other_node = p_vessel->GetNodeAtOppositeEnd(p_node);
            typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter = std::find(nodes.begin(), nodes.end(), p_other_node);
            unsigned other_node_index = std::distance(nodes.begin(), node_iter);
            node_indexes.push_back(other_node_index);
        }
        connectivity.push_back(node_indexes);
    }
    return connectivity;
}

template <unsigned DIM>
std::vector<std::vector<unsigned> > VascularNetwork<DIM>::GetNodeVesselConnectivity()
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels = GetVessels();
    unsigned num_nodes = nodes.size();
    std::vector<std::vector<unsigned> > connectivity;
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        boost::shared_ptr<VascularNode<DIM> > p_node = nodes[node_index];
        std::vector<unsigned> vessel_indexes;

        unsigned num_segments_on_node = p_node->GetNumberOfSegments();
        for (unsigned segment_index = 0; segment_index < num_segments_on_node; segment_index++)
        {
            boost::shared_ptr<Vessel<DIM> > p_vessel = p_node->GetVesselSegment(segment_index)->GetVessel();

            typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator vessel_iter =
                    std::find(vessels.begin(), vessels.end(), p_vessel);
            unsigned vessel_index = std::distance(vessels.begin(), vessel_iter);
            vessel_indexes.push_back(vessel_index);
        }
        connectivity.push_back(vessel_indexes);
    }
    return connectivity;
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetMaxBranchesOnNode()
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();
    unsigned num_nodes = nodes.size();

    // Get maximum number of segments attached to a node in the whole network.
    unsigned max_num_branches = 0;
    for(unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        boost::shared_ptr<VascularNode<DIM> > p_each_node = nodes[node_index];
        unsigned num_segments_on_node = nodes[node_index]->GetNumberOfSegments();

        if (num_segments_on_node > max_num_branches)
        {
            max_num_branches = num_segments_on_node;
        }
    }
    return max_num_branches;
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetVesselIndex(boost::shared_ptr<Vessel<DIM> > pVessel)
{
    unsigned index = 0;
    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        if(*it == pVessel)
        {
            return index;
        }
        index ++;
    }
    EXCEPTION("Input vessel is not in the network.");
}

template <unsigned DIM>
unsigned VascularNetwork<DIM>::GetVesselSegmentIndex(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment)
{
    unsigned index = 0;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        if(*it == pVesselSegment)
        {
            return index;
        }
        index++;
    }
    EXCEPTION("Input vessel is not in the network.");
}

template <unsigned DIM>
std::vector<boost::shared_ptr<VesselSegment<DIM> > > VascularNetwork<DIM>::GetVesselSegments()
{
    if(!mSegmentsUpToDate)
    {
        UpdateSegments();
    }

    return mSegments;
}

/*
 * Helper class for "connected" methods
 */
template<typename TimeMap> class bfs_time_visitor : public boost::default_bfs_visitor
{
    typedef typename boost::property_traits<TimeMap>::value_type T;

public:

    TimeMap m_timemap;
    T& m_time;

    bfs_time_visitor(TimeMap tmap, T& t)
    :m_timemap(tmap),
     m_time(t)
    {
    }

    template<typename Vertex, typename Graph>
    void discover_vertex(Vertex u, const Graph& g) const
    {
        put(m_timemap, u, m_time++);
    }
};

template <unsigned DIM>
bool VascularNetwork<DIM>::IsConnected(boost::shared_ptr<VascularNode<DIM> > pSourceNode, boost::shared_ptr<VascularNode<DIM> > pQueryNode)
{

    if (!NodeIsInNetwork(pSourceNode))
    {
        EXCEPTION("Source node is not in network.");
    }
    if (!NodeIsInNetwork(pQueryNode))
    {
        EXCEPTION("Query node is not in network.");
    }

    if (pSourceNode == pQueryNode)
    {
        return true;
    }

    boost::shared_ptr<Vessel<DIM> > p_source_vessel = pSourceNode->GetVesselSegment(0)->GetVessel();
    boost::shared_ptr<Vessel<DIM> > p_query_vessel = pQueryNode->GetVesselSegment(0)->GetVessel();

    if (p_source_vessel == p_query_vessel || p_source_vessel->IsConnectedTo(p_query_vessel))
    {
        return true;
    }

    // Assign the vessel nodes unique IDs
    std::vector<boost::shared_ptr<VascularNode<DIM> > >  vessel_nodes = GetVesselEndNodes();
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    unsigned counter = 0;
    for(node_iter = vessel_nodes.begin(); node_iter != vessel_nodes.end(); node_iter++)
    {
        (*node_iter)->SetId(counter);
        counter ++;
    }

    // Get the start nodes of the containing and query vessels
    boost::shared_ptr<VascularNode<DIM> > pEquivalentSourceNode = p_source_vessel->GetStartNode();
    boost::shared_ptr<VascularNode<DIM> > pEquivalentQueryNode = p_query_vessel->GetStartNode();

    // construct graph representation of vessel network
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < mVessels.size(); i++)
    {
        add_edge(mVessels[i]->GetStartNode()->GetId(), mVessels[i]->GetEndNode()->GetId(), G);
    }

    // typedefs
    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    // a vector to hold the discover time property for each vertex
    std::vector<Size> dtime(num_vertices(G));

    Size time = 0;
    bfs_time_visitor<Size*>vis(&dtime[0], time);

    // use breadth first search to establish discovery time of all nodes from node1
    // this assigns a discovery time to dTime for each node (index relates to nodeID)
    // dTime is zero for node1 and all other nodes that are not connected to node1
    // dTime is nonzero for all nodes that are connected to node1 (except node 1 itself)
    breadth_first_search(G,vertex(pEquivalentSourceNode->GetId(),G), boost::visitor(vis));

    return (dtime[pEquivalentQueryNode->GetId()] > 0);
}

template <unsigned DIM>
std::vector<bool > VascularNetwork<DIM>::IsConnected(std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes,
        std::vector<boost::shared_ptr<VascularNode<DIM> > > queryNodes)
        {
    // Assign the vessel nodes unique IDs
    std::vector<boost::shared_ptr<VascularNode<DIM> > >  vessel_nodes = GetVesselEndNodes();
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    unsigned counter = 0;
    for(node_iter = vessel_nodes.begin(); node_iter != vessel_nodes.end(); node_iter++)
    {
        (*node_iter)->SetId(counter);
        counter ++;
    }

    // construct graph representation of vessel network
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    Graph G;

    for (unsigned i = 0; i < mVessels.size(); i++)
    {
        add_edge(mVessels[i]->GetStartNode()->GetId(), mVessels[i]->GetEndNode()->GetId(), G);
    }

    std::vector<bool > connected(queryNodes.size(),false);

    for(unsigned i=0; i<sourceNodes.size(); i++)
    {

        if (!NodeIsInNetwork(sourceNodes[i]))
        {
            EXCEPTION("Source node is not in network.");
        }

        boost::shared_ptr<VascularNode<DIM> > pSourceNode = sourceNodes[i];
        boost::shared_ptr<Vessel<DIM> > p_source_vessel = pSourceNode->GetVesselSegment(0)->GetVessel();
        boost::shared_ptr<VascularNode<DIM> > pEquivalentSourceNode = p_source_vessel->GetStartNode();

        // a vector to hold the discover time property for each vertex
        std::vector<Size> dtime(num_vertices(G));

        Size time = 0;
        bfs_time_visitor<Size*>vis(&dtime[0], time);

        breadth_first_search(G,vertex(pEquivalentSourceNode->GetId(),G), boost::visitor(vis));

        for (unsigned j=0; j<queryNodes.size(); j++)
        {

            if (!NodeIsInNetwork(queryNodes[j]))
            {
                EXCEPTION("Query node is not in network.");
            }

            boost::shared_ptr<VascularNode<DIM> > pQueryNode = queryNodes[j];

            if (pSourceNode == pQueryNode)
            {
                connected[j] = true;
                continue;
            }

            boost::shared_ptr<Vessel<DIM> > p_query_vessel = pQueryNode->GetVesselSegment(0)->GetVessel();
            if (p_source_vessel == p_query_vessel || p_source_vessel->IsConnectedTo(p_query_vessel))
            {
                connected[j] = true;
                continue;
            }

            boost::shared_ptr<VascularNode<DIM> > pEquivalentQueryNode = p_query_vessel->GetStartNode();

            if (dtime[pEquivalentQueryNode->GetId()] > 0)
            {
                connected[j] = true;
            }
        }

    }
    return connected;
        }

template <unsigned DIM>
bool VascularNetwork<DIM>::NodeIsInNetwork(boost::shared_ptr<VascularNode<DIM> > pSourceNode)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    return (std::find(nodes.begin(), nodes.end(), pSourceNode) != nodes.end());
}

template <unsigned DIM>
void VascularNetwork<DIM>::MergeCoincidentNodes(double tolerance)
{
    MergeCoincidentNodes(GetNodes(), tolerance);
}

template <unsigned DIM>
void VascularNetwork<DIM>::MergeCoincidentNodes(std::vector<boost::shared_ptr<Vessel<DIM> > > pVessels, double tolerance)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    for(unsigned idx = 0; idx <pVessels.size(); idx++)
    {
        std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = pVessels[idx]->GetNodes();
        nodes.insert(nodes.end(), vessel_nodes.begin(), vessel_nodes.end());
    }
    MergeCoincidentNodes(nodes, tolerance);
}

template <unsigned DIM>
void VascularNetwork<DIM>::MergeCoincidentNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes, double tolerance)
{
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it;
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it2;
    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it3;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        for(it2 = nodes.begin(); it2 != nodes.end(); it2++)
        {
            // If the nodes are not identical
            if ((*it) != (*it2))
            {
                // If the node locations are the same - according to the ChastePoint definition
                bool is_coincident = false;
                if(tolerance >0.0)
                {
                    is_coincident = (*it)->GetDistance((*it2)) <= tolerance;
                }
                else
                {
                    is_coincident = (*it)->IsCoincident((*it2));
                }

                if(is_coincident)
                {
                    // Replace the node corresponding to 'it2' with the one corresponding to 'it'
                    // in all segments.
                    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = (*it2)->GetVesselSegments();
                    for(it3 = segments.begin(); it3 != segments.end(); it3++)
                    {
                        if ((*it3)->GetNode(0) == (*it2))
                        {
                            (*it3)->ReplaceNode(0, (*it));
                        }
                        else if(((*it3)->GetNode(1) == (*it2)))
                        {
                            (*it3)->ReplaceNode(1, (*it));
                        }
                    }
                }
            }
        }
    }
    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void VascularNetwork<DIM>::SetNodeData(VasculatureData data)
{
    //NEVER_REACHED;
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::SetVesselData(VasculatureData data)
{
    //NEVER_REACHED;
    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::SetSegmentData(VasculatureData data)
{

    //NEVER_REACHED;
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::SetSegmentProperties(boost::shared_ptr<VesselSegment<DIM> >  prototype)
{

    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    typename std::vector<boost::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetRadius(prototype->GetRadius());
        (*it)->GetFlowProperties()->SetImpedance(prototype->GetFlowProperties()->GetImpedance());
        (*it)->GetFlowProperties()->SetHaematocrit(prototype->GetFlowProperties()->GetHaematocrit());
        (*it)->GetFlowProperties()->SetFlowRate(prototype->GetFlowProperties()->GetFlowRate());
        (*it)->GetFlowProperties()->SetViscosity(prototype->GetFlowProperties()->GetViscosity());
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::Translate(const c_vector<double, DIM>& rTranslationVector)
{
    Translate(rTranslationVector, mVessels);
}

template <unsigned DIM>
void VascularNetwork<DIM>::Translate(const c_vector<double, DIM>& rTranslationVector, std::vector<boost::shared_ptr<Vessel<DIM> > > vessels)
{
    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes;
    for(unsigned idx = 0; idx <vessels.size(); idx++)
    {
        std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = vessels[idx]->GetNodes();
        std::copy(vessel_nodes.begin(), vessel_nodes.end(), std::inserter(nodes, nodes.begin()));
    }

    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {
        (*node_iter)->SetLocation((*node_iter)->GetLocationVector() + rTranslationVector);
    }
}

template <unsigned DIM>
void VascularNetwork<DIM>::RemoveVessel(boost::shared_ptr<Vessel<DIM> > pVessel, bool deleteVessel)
{

    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it = std::find(mVessels.begin(), mVessels.end(), pVessel);
    if(it != mVessels.end())
    {
        if(deleteVessel)
        {
            (*it)->Remove();
            // todo remove cells from nodes if the nodes will be removed from the network
        }
        mVessels.erase(it);
    }
    else
    {
        EXCEPTION("Vessel is not contained inside network.");
    }

    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void VascularNetwork<DIM>::SetNodeRadii(double radius)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();

    for(unsigned idx=0; idx<nodes.size();idx++)
    {
        nodes[idx]->SetRadius(radius);
    }

}

template <unsigned DIM>
void VascularNetwork<DIM>::SetSegmentRadii(double radius)
{
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = GetVesselSegments();

    for(unsigned idx=0; idx<segments.size();idx++)
    {
        segments[idx]->SetRadius(radius);
    }
}

#ifdef CHASTE_VTK
template <unsigned DIM>
vtkSmartPointer<vtkPolyData> VascularNetwork<DIM>::GetVtk()
{
    UpdateVesselIds();

    // Set up the vessel and node data arrays.
    std::vector<vtkSmartPointer<vtkDoubleArray> > pVesselInfoVector;
    std::map<std::string, boost::any >::iterator map_iterator;
    std::map<std::string, boost::any > generic_data_map = mVessels[0]->rGetDataContainer().GetMap();
    for(map_iterator = generic_data_map.begin(); map_iterator != generic_data_map.end(); map_iterator++)
    {
        vtkSmartPointer<vtkDoubleArray> pVesselInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pVesselInfo->SetNumberOfComponents(1);
        pVesselInfo->SetNumberOfTuples(mVessels.size());
        pVesselInfo->SetName((*map_iterator).first.c_str());

        // Only write information that can be cast to double
        if (map_iterator->second.type() == typeid(double) || map_iterator->second.type() == typeid(unsigned) ||
                map_iterator->second.type() == typeid(bool))
        {
            pVesselInfoVector.push_back(pVesselInfo);
        }
    }
    std::map<std::string, double>::iterator vessel_map_iterator;
    std::map<std::string, double> vessel_data_map = mVessels[0]->GetVtkData();
    for(vessel_map_iterator = vessel_data_map.begin(); vessel_map_iterator != vessel_data_map.end(); vessel_map_iterator++)
    {
        vtkSmartPointer<vtkDoubleArray> pVesselInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pVesselInfo->SetNumberOfComponents(1);
        pVesselInfo->SetNumberOfTuples(mVessels.size());
        pVesselInfo->SetName((*vessel_map_iterator).first.c_str());
        pVesselInfoVector.push_back(pVesselInfo);
    }

    // Set up node info arrays
    std::vector<vtkSmartPointer<vtkDoubleArray> > pNodeInfoVector;
    unsigned numberOfNodes = this->GetNumberOfNodes();

    if (mVessels.size()==0)
    {
        EXCEPTION("Attempting to write a network with no vessels");
    }

    std::map<std::string, boost::any>::iterator generic_node_map_iterator;
    std::map<std::string, boost::any> generic_node_map = mVessels[0]->GetStartNode()->rGetDataContainer().GetMap();
    for(generic_node_map_iterator = generic_node_map.begin(); generic_node_map_iterator != generic_node_map.end(); generic_node_map_iterator++)
    {
        vtkSmartPointer<vtkDoubleArray> pNodeInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pNodeInfo->SetNumberOfComponents(1);
        pNodeInfo->SetNumberOfTuples(numberOfNodes);
        pNodeInfo->SetName((*generic_node_map_iterator).first.c_str());

        // Only write information that can be cast to double
        if (generic_node_map_iterator->second.type() == typeid(double) || generic_node_map_iterator->second.type() == typeid(unsigned) ||
                generic_node_map_iterator->second.type() == typeid(bool))
        {
            pNodeInfoVector.push_back(pNodeInfo);
        }
    }

    std::map<std::string, double>::iterator vtk_node_map_iterator;
    std::map<std::string, double> vtk_node_map = mVessels[0]->GetStartNode()->GetVtkData();
    for(vtk_node_map_iterator = vtk_node_map.begin(); vtk_node_map_iterator != vtk_node_map.end(); vtk_node_map_iterator++)
    {
        vtkSmartPointer<vtkDoubleArray> pNodeInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pNodeInfo->SetNumberOfComponents(1);
        pNodeInfo->SetNumberOfTuples(numberOfNodes);
        pNodeInfo->SetName((*vtk_node_map_iterator).first.c_str());
        pNodeInfoVector.push_back(pNodeInfo);
    }

    // Create the geometric data
    vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pPoints= vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> pLines = vtkSmartPointer<vtkCellArray>::New();

    unsigned vessel_index=0;

    // Do the nodes
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = this->GetNodes();
    for(unsigned idx=0; idx<nodes.size(); idx++)
    {
        nodes[idx]->SetTempId(idx);
        if(DIM == 2)
        {
            pPoints->InsertNextPoint(nodes[idx]->GetLocationVector()[0], nodes[idx]->GetLocationVector()[1], 0.0);
        }
        else
        {
            pPoints->InsertNextPoint(nodes[idx]->GetLocationVector()[0], nodes[idx]->GetLocationVector()[1], nodes[idx]->GetLocationVector()[2]);
        }

        std::map<std::string, double> vtk_node_data = nodes[idx]->GetVtkData();
        std::map<std::string, boost::any> generic_node_data = nodes[idx]->rGetDataContainer().GetMap();

        // Add the node data
        for(unsigned jdx=0; jdx < pNodeInfoVector.size(); jdx++)
        {
            // Get the key
            std::string key = pNodeInfoVector[jdx]->GetName();

            // If it is in the vtk data use it
            if(vtk_node_data.count(key) == 1)
            {
                pNodeInfoVector[jdx]->SetValue(idx, vtk_node_data[key]);
            }
            // Otherwise check the generic data
            else if(generic_node_data.count(key) == 1)
            {
                if(generic_node_data[key].type() == typeid(double))
                {
                    double cast_value = boost::any_cast<double>(generic_node_data[key]);
                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
                }
                else if(generic_node_data[key].type() == typeid(unsigned))
                {
                    double cast_value = double(boost::any_cast<unsigned>(generic_node_data[key]));
                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
                }
                else if(generic_node_data[key].type() == typeid(bool))
                {
                    double cast_value = double(boost::any_cast<bool>(generic_node_data[key]));
                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
                }
            }
        }
    }

    // Do the vessels
    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it < mVessels.end(); it++)
    {
        vtkSmartPointer<vtkLine> pLine = vtkSmartPointer<vtkLine>::New();
        std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = (*it)->GetSegments();

        for(unsigned i = 0; i < segments.size(); i++)
        {
            pLine->GetPointIds()->InsertId(i, segments[i]->GetNode(0)->GetTempId());

            // Do an extra insert for the last node in the segment
            if (i == segments.size() - 1)
            {
                pLine->GetPointIds()->InsertId(i + 1, segments[i]->GetNode(1)->GetTempId());
            }
        }
        pLines->InsertNextCell(pLine);

        // Add the vessel data
        std::map<std::string, double> vtk_vessel_data = (*it)->GetVtkData();
        std::map<std::string, boost::any> generic_vessel_data = (*it)->rGetDataContainer().GetMap();

        for(unsigned idx=0; idx < pVesselInfoVector.size(); idx++)
        {
            // Get the key
            std::string key = pVesselInfoVector[idx]->GetName();
            // If it is in the vtk data use it
            if(vtk_vessel_data.count(key) == 1)
            {
                pVesselInfoVector[idx]->SetValue(vessel_index, vtk_vessel_data[key]);
            }
            // Otherwise check the generic data
            else if(generic_vessel_data.count(key) == 1)
            {
                if(generic_vessel_data[key].type() == typeid(double))
                {
                    double cast_value = boost::any_cast<double>(generic_vessel_data[key]);
                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
                }
                else if(generic_vessel_data[key].type() == typeid(unsigned))
                {
                    double cast_value = double(boost::any_cast<unsigned>(generic_vessel_data[key]));
                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
                }
                else if(generic_vessel_data[key].type() == typeid(bool))
                {
                    double cast_value = double(boost::any_cast<bool>(generic_vessel_data[key]));
                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
                }
            }
        }
        vessel_index++;
    }
    pPolyData->SetPoints(pPoints);
    pPolyData->SetLines(pLines);

    for (unsigned i = 0; i < pVesselInfoVector.size(); i++)
    {
        pPolyData->GetCellData()->AddArray(pVesselInfoVector[i]);
    }

    for (unsigned i = 0; i < pNodeInfoVector.size(); i++)
    {
        pPolyData->GetPointData()->AddArray(pNodeInfoVector[i]);
    }
    return pPolyData;
}

template <unsigned DIM>
void VascularNetwork<DIM>::Write(const std::string& filename)
{
    vtkSmartPointer<vtkPolyData> p_polydata = GetVtk();
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(p_polydata);
    writer->Write();
}
#endif // CHASTE_VTK

template <unsigned DIM>
void VascularNetwork<DIM>::WriteConnectivity(const std::string& output_filename)
{
    // construct graph representation of vessel network
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iterator;
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();

    for (node_iterator = nodes.begin(); node_iterator != nodes.end(); node_iterator++)
    {
        if ((*node_iterator)->GetVesselSegments().size() > 1)
        {
            for (unsigned j = 1; j < (*node_iterator)->GetVesselSegments().size(); j++)
            {
                add_edge(GetVesselIndex((*node_iterator)->GetVesselSegment(0)->GetVessel()),
                        GetVesselIndex((*node_iterator)->GetVesselSegment(j)->GetVessel()), G);
            }
        }
    }

    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator vessel_iterator;
    for (vessel_iterator = mVessels.begin(); vessel_iterator != mVessels.end(); vessel_iterator++)
    {
        if ((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
        {
            add_edge(GetVesselIndex((*vessel_iterator)),
                    GetVesselIndex((*vessel_iterator)), G);
        }
    }

    std::ofstream outf(output_filename.c_str());
    boost::dynamic_properties dp;
    dp.property("node_id", get(boost::vertex_index, G));
    write_graphviz_dp(outf, G, dp);
}

template <unsigned DIM>
void VascularNetwork<DIM>::UpdateAll(bool merge)
{
    if(merge)
    {
        MergeCoincidentNodes();
    }
    UpdateSegments();
    UpdateVesselNodes();
    UpdateNodes();
    UpdateVesselIds();
}

template<unsigned DIM>
struct NodePtrComp
{
  bool operator()( const boost::shared_ptr<VascularNode<DIM> >  & a, const boost::shared_ptr<VascularNode<DIM> >  & b )
    { return a->GetTempId() > b->GetTempId(); }
};

template<unsigned DIM>
void VascularNetwork<DIM>::UpdateNodes()
{
      mNodes = std::vector<boost::shared_ptr<VascularNode<DIM> > >();
      std::set<boost::shared_ptr<VascularNode<DIM> >, NodePtrComp<DIM> >  nodes;
      std::vector<boost::shared_ptr<VascularNode<DIM> > > temp_nodes;

      typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;

      unsigned counter = 0;
      for(it = mVessels.begin(); it != mVessels.end(); it++)
      {
          std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = (*it)->GetNodes();
          for (unsigned idx=0; idx<vessel_nodes.size(); idx++)
          {
              vessel_nodes[idx]->SetTempId(counter);
              temp_nodes.push_back(vessel_nodes[idx]);
              counter ++;
          }
      }
      std::copy(temp_nodes.begin(), temp_nodes.end(), std::inserter(nodes, nodes.begin()));

      std::copy(nodes.begin(), nodes.end(), std::back_inserter(mNodes));
      mNodesUpToDate = true;


//    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
//
//    for(it = mVessels.begin(); it != mVessels.end(); it++)
//    {
//        std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = (*it)->GetNodes();
//        std::copy(vessel_nodes.begin(), vessel_nodes.end(), std::inserter(nodes, nodes.begin()));
//    }
//    std::copy(nodes.begin(), nodes.end(), std::back_inserter(mNodes));
//    mNodesUpToDate = true;
}

template<unsigned DIM>
void VascularNetwork<DIM>::UpdateSegments()
{
    mSegments = std::vector<boost::shared_ptr<VesselSegment<DIM> > >();
    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        std::vector<boost::shared_ptr<VesselSegment<DIM> > > vessel_segments = (*it)->GetSegments();
        std::copy(vessel_segments.begin(), vessel_segments.end(), std::back_inserter(mSegments));
    }
    mSegmentsUpToDate = true;
}

template<unsigned DIM>
void VascularNetwork<DIM>::UpdateVesselNodes()
{
    mVesselNodes = std::vector<boost::shared_ptr<VascularNode<DIM> > >();
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        (*it)->UpdateNodes();
        boost::shared_ptr<VascularNode<DIM> > vessel_node1 = (*it)->GetStartNode();
        nodes.insert(vessel_node1);
        boost::shared_ptr<VascularNode<DIM> > vessel_node2 = (*it)->GetEndNode();
        nodes.insert(vessel_node2);
    }
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(mVesselNodes));
    mVesselNodesUpToDate = true;
}

template<unsigned DIM>
void VascularNetwork<DIM>::UpdateVesselIds()
{
    for(unsigned idx=0;idx<mVessels.size();idx++)
    {
        mVessels[idx]->SetId(idx);
    }
}

template<unsigned DIM>
bool VascularNetwork<DIM>::VesselCrossesLineSegment(c_vector<double, DIM> coordinate_1, c_vector<double, DIM> coordinate_2, double radius)
{

    boost::shared_ptr<VesselSegment<DIM> > temp_segment = VesselSegment<DIM>::Create(VascularNode<DIM>::Create(coordinate_1), VascularNode<DIM>::Create(coordinate_2));

    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> nearest_segment = GetNearestSegment(temp_segment);

    // todo a false here does not necessarily guarantee that a vessel does not cross a line segment since get nearest
    // segment only returns one segment
    bool crosses_segment = (nearest_segment.second <= radius) && (nearest_segment.first->GetDistance(coordinate_1) > radius) && (nearest_segment.first->GetDistance(coordinate_2) > radius);

    return  crosses_segment;

}

// Explicit instantiation
template class VascularNetwork<2>;
template class VascularNetwork<3>;

