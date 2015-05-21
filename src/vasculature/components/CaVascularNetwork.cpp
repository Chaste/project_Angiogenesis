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

#include "CaVascularNetwork.hpp"

template <unsigned DIM>
CaVascularNetwork<DIM>::CaVascularNetwork()
: mVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > >()),
  mSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >()),
  mSegmentsUpToDate(false),
  mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
  mNodesUpToDate(false),
  mVesselNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
  mVesselNodesUpToDate(false),
  mDataContainer()
  {
  }

template <unsigned DIM>
CaVascularNetwork<DIM>::~CaVascularNetwork()
{
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel)
{
    mVessels.push_back(pVessel);
    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels)
{
    mVessels.insert(mVessels.end(), vessels.begin(), vessels.end());
    mSegmentsUpToDate = false;
    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::CopyVessels()
{
    return CopyVessels(mVessels);
}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::CopyVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels)
{
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iter;
    std::vector<boost::shared_ptr<CaVessel<DIM> > > new_vessels;
    for(vessel_iter = vessels.begin(); vessel_iter != vessels.end(); vessel_iter++)
    {
        typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter;
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > new_segments;
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*vessel_iter)->GetSegments();
        for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
        {
            ChastePoint<DIM> node0_location = (*segment_iter)->GetNode(0)->GetLocation();
            ChastePoint<DIM> node1_location = (*segment_iter)->GetNode(1)->GetLocation();
            new_segments.push_back(CaVesselSegment<DIM>::Create(VascularNode<DIM>::Create(node0_location),
                                                                VascularNode<DIM>::Create(node1_location)));
        }
        new_vessels.push_back(CaVessel<DIM>::Create(new_segments));
    }

    MergeCoincidentNodes(new_vessels);
    AddVessels(new_vessels);
    return new_vessels;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::RemoveVessel(boost::shared_ptr<CaVessel<DIM> > pVessel)
{

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it = std::find(mVessels.begin(), mVessels.end(), pVessel);
    if(it != mVessels.end())
    {
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
std::vector<std::pair<double, double> > CaVascularNetwork<DIM>::GetExtents()
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
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::GetNearestNode(const ChastePoint<DIM>& rLocation)
{
    return GetNearestNode(rLocation.rGetLocation());
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::GetNearestNode(c_vector<double, DIM> location)
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
std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double>  CaVascularNetwork<DIM>::GetNearestSegment(const ChastePoint<DIM>& rLocation)
{
    return GetNearestSegment(rLocation.rGetLocation());
}

template <unsigned DIM>
std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double>  CaVascularNetwork<DIM>::GetNearestSegment(c_vector<double, DIM> location)
{
    boost::shared_ptr<CaVesselSegment<DIM> > nearest_segment;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = GetVesselSegments();

    double min_distance = 1.e12;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter;
    for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
    {
        double segment_distance = (*segment_iter)->GetDistance(location);
        if (segment_distance < min_distance)
        {
            min_distance = segment_distance;
            nearest_segment = (*segment_iter) ;
        }
    }
    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> return_pair =
            std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double>(nearest_segment, min_distance);
    return return_pair;
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetNearestVessel(const ChastePoint<DIM>& rLocation)
{
    return GetNearestSegment(rLocation).first->GetVessel();
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetNearestVessel(c_vector<double, DIM> location)
{
    return GetNearestSegment(location).first->GetVessel();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::NumberOfNodesNearLocation(const ChastePoint<DIM>& rLocation, double radius)
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
std::vector<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::GetNodes()
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::GetNode(unsigned index)
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes[index];
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfNodes()
{
    if(!mNodesUpToDate)
    {
        UpdateNodes();
    }

    return mNodes.size();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfVesselNodes()
{
    if(!mVesselNodesUpToDate)
    {
        UpdateVesselNodes();
    }

    return mVesselNodes.size();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNodeIndex(boost::shared_ptr<VascularNode<DIM> > node)
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
unsigned CaVascularNetwork<DIM>::GetNumberOfVessels()
{
    return mVessels.size();
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetVessel(unsigned index)
{
    return mVessels[index];
}

template <unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::GetVesselEndNodes()
{

    if(!mVesselNodesUpToDate)
    {
        UpdateVesselNodes();
    }

    return mVesselNodes;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::GetVessels()
{
    return mVessels;
}

template <unsigned DIM>
std::vector<std::vector<unsigned> > CaVascularNetwork<DIM>::GetNodeNodeConnectivity()
{
	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();
	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = GetVessels();
	std::vector<std::vector<unsigned> > node_vessel_connectivity = GetNodeVesselConnectivity();

	std::vector<std::vector<unsigned> > connectivity;
	for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
	{
		std::vector<unsigned> node_indexes;
		boost::shared_ptr<VascularNode<DIM> > p_node = nodes[node_index];
		unsigned num_branches = node_vessel_connectivity[node_index].size();
		for (unsigned vessel_index = 0; vessel_index < num_branches; vessel_index++)
		{
			boost::shared_ptr<CaVessel<DIM> > p_vessel = vessels[node_vessel_connectivity[node_index][vessel_index]];

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
std::vector<std::vector<unsigned> > CaVascularNetwork<DIM>::GetNodeVesselConnectivity()
{
	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();
	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = GetVessels();
	unsigned num_nodes = nodes.size();
	std::vector<std::vector<unsigned> > connectivity;

	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		boost::shared_ptr<VascularNode<DIM> > p_node = nodes[node_index];
		std::vector<unsigned> vessel_indexes;

		unsigned num_segments_on_node = p_node->GetNumberOfSegments();
		for (unsigned segment_index = 0; segment_index < num_segments_on_node; segment_index++)
		{
			boost::shared_ptr<CaVessel<DIM> > p_vessel = p_node->GetVesselSegment(segment_index)->GetVessel();

			typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iter =
					std::find(vessels.begin(), vessels.end(), p_vessel);
			unsigned vessel_index = std::distance(vessels.begin(), vessel_iter);
			vessel_indexes.push_back(vessel_index);
		}
		connectivity.push_back(vessel_indexes);
	}
	return connectivity;
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetMaxBranchesOnNode()
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
unsigned CaVascularNetwork<DIM>::GetVesselIndex(boost::shared_ptr<CaVessel<DIM> > pVessel)
{
    unsigned index = 0;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
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
unsigned CaVascularNetwork<DIM>::GetVesselSegmentIndex(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment)
{
    unsigned index = 0;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = GetVesselSegments();
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it;
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
std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > CaVascularNetwork<DIM>::GetVesselSegments()
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
bool CaVascularNetwork<DIM>::IsConnected(boost::shared_ptr<VascularNode<DIM> > pSourceNode, boost::shared_ptr<VascularNode<DIM> > pQueryNode)
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

    boost::shared_ptr<CaVessel<DIM> > p_source_vessel = pSourceNode->GetVesselSegment(0)->GetVessel();
    boost::shared_ptr<CaVessel<DIM> > p_query_vessel = pQueryNode->GetVesselSegment(0)->GetVessel();

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
std::vector<bool > CaVascularNetwork<DIM>::IsConnected(std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes,
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
        boost::shared_ptr<CaVessel<DIM> > p_source_vessel = pSourceNode->GetVesselSegment(0)->GetVessel();
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

            boost::shared_ptr<CaVessel<DIM> > p_query_vessel = pQueryNode->GetVesselSegment(0)->GetVessel();
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
bool CaVascularNetwork<DIM>::NodeIsInNetwork(boost::shared_ptr<VascularNode<DIM> > pSourceNode)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    return (std::find(nodes.begin(), nodes.end(), pSourceNode) != nodes.end());
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::MergeCoincidentNodes()
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    MergeCoincidentNodes(nodes);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::MergeCoincidentNodes(std::vector<boost::shared_ptr<CaVessel<DIM> > > pVessels)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    for(unsigned idx = 0; idx <pVessels.size(); idx++)
    {
        std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = pVessels[idx]->GetNodes();
        nodes.insert(nodes.end(), vessel_nodes.begin(), vessel_nodes.end());
    }
    MergeCoincidentNodes(nodes);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::MergeCoincidentNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes)
{
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it;
    typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator it2;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it3;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        for(it2 = nodes.begin(); it2 != nodes.end(); it2++)
        {
            // If the nodes are not identical
            if ((*it) != (*it2))
            {
                // If the node locations are the same - according to the ChastePoint definition
                if((*it)->IsCoincident((*it2)))
                {
                    // Replace the node corresponding to 'it2' with the one corresponding to 'it'
                    // in all segments.
                    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it2)->GetVesselSegments();
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

    mNodesUpToDate = false;
    mVesselNodesUpToDate = false;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetNodeData(VasculatureData data)
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
void CaVascularNetwork<DIM>::SetVesselData(VasculatureData data)
{
    //NEVER_REACHED;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetSegmentData(VasculatureData data)
{

    //NEVER_REACHED;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = GetVesselSegments();

    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetSegmentProperties(boost::shared_ptr<CaVesselSegment<DIM> >  prototype)
{

    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = GetVesselSegments();

    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it;
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
void CaVascularNetwork<DIM>::Translate(const c_vector<double, DIM>& rTranslationVector)
{
    Translate(rTranslationVector, mVessels);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::Translate(const c_vector<double, DIM>& rTranslationVector, std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels)
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
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::DivideVessel(boost::shared_ptr<CaVessel<DIM> > pVessel, ChastePoint<DIM> location)
{

    // assert that location must not soley be on the end of a vessel
    // it may, however, be on the end of a vessel and lie in the middle of the vessel too

    boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment;

    if (pVessel->GetStartNode()->IsCoincident(location) || pVessel->GetEndNode()->IsCoincident(location) )
    {
        unsigned locatedInsideVessel = 0;

        for (unsigned i = 0; i < pVessel->GetNumberOfSegments(); i++)
        {
            if (pVessel->GetSegment(i)->GetDistance(location) <= 1e-6)
            {
                locatedInsideVessel++;
                if (i != 0 && i != pVessel->GetNumberOfSegments() - 1)
                {

                    pVesselSegment = pVessel->GetSegment(i);
                }
            }
        }

        // vessel may only ever be allowed to exist twice in the same location
        // this is also only applicable when a vessel is looping around on itself
        if(locatedInsideVessel != 2)
        {
            EXCEPTION("If the location of a divide coincides with the end of a vessel then there must be at least one other segment in the centre of that vessel at which the vessel can be divided. A vessel cannot be divided at its end.");
        }
        // both vessel nodes cannot be present at the same location in this case
        if(pVessel->GetStartNode()->IsCoincident(pVessel->GetEndNode()))
        {
            EXCEPTION("If the location of a divide is on the end of a vessel the two ends of a vessel cannot coincide.");
        }
    }
    else
    {
        bool locatedInsideVessel = false;

        for (unsigned i = 0; i < pVessel->GetNumberOfSegments(); i++)
        {
            if (pVessel->GetSegment(i)->GetDistance(location) <= 1e-6)
            {
                locatedInsideVessel = true;
                pVesselSegment = pVessel->GetSegment(i);
                break;
            }
        }

        if(!locatedInsideVessel)
        {
            EXCEPTION("The division cannot take place at the location specified since the vessel does not exist at that location.");
        }
    }

    assert(pVesselSegment);

    boost::shared_ptr<VascularNode<DIM> > p_new_node = pVessel->DivideSegment(location);

    // create two new vessels and give them the properties of the vessel to be divided and appropriately divide the vessel segments of the divided vessel between the two new vessels
    // new vessels will share a new common vessel network node
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > start_segments;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > end_segments;
    unsigned segment_index;
    for (unsigned i = 0; i < pVessel->GetNumberOfSegments(); i++)
    {
        start_segments.push_back(pVessel->GetSegment(i));
        if (pVessel->GetSegment(i)->GetNode(0)->IsCoincident(location) )
        {
            segment_index = i;
            break;
        }
    }

    for (unsigned i = segment_index; i < pVessel->GetNumberOfSegments(); i++)
    {
        end_segments.push_back(pVessel->GetSegment(i));
    }
    boost::shared_ptr<CaVessel<DIM> > p_new_vessel1 = CaVessel<DIM>::Create(start_segments);
    boost::shared_ptr<CaVessel<DIM> > p_new_vessel2 = CaVessel<DIM>::Create(end_segments);
    p_new_vessel1->CopyDataFromExistingVessel(pVessel);
    p_new_vessel2->CopyDataFromExistingVessel(pVessel);

    AddVessel(p_new_vessel1);
    AddVessel(p_new_vessel2);
    RemoveVessel(pVessel);

    return p_new_node;
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::FormSprout(ChastePoint<DIM> sproutBaseLocation, ChastePoint<DIM> sproutTipLocation)
{

    // locate vessel at which the location of the sprout base exists

    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> nearest_segment = GetNearestSegment(sproutBaseLocation);
    if (nearest_segment.second > 1e-6)
    {
        EXCEPTION("No vessel located at sprout base.");
    }

    // divide vessel at location of sprout base

    boost::shared_ptr<VascularNode<DIM> > p_new_node = DivideVessel(nearest_segment.first->GetVessel(), sproutBaseLocation);

    if(p_new_node->IsMigrating())
    {
        EXCEPTION("Cannot form sprout at a location which is migrating (we assume migrating nodes are occupied by a tip cell).");
    }

    // create new vessel

    boost::shared_ptr<VascularNode<DIM> > p_new_node_at_tip = VascularNode<DIM>::Create(p_new_node);
    p_new_node_at_tip->SetLocation(sproutTipLocation);
    p_new_node_at_tip->SetIsMigrating(true);
    boost::shared_ptr<CaVesselSegment<DIM> > p_new_segment = CaVesselSegment<DIM>::Create(p_new_node,p_new_node_at_tip);
    p_new_segment->CopyDataFromExistingSegment(nearest_segment.first);
    boost::shared_ptr<CaVessel<DIM> > p_new_vessel = CaVessel<DIM>::Create(p_new_segment);
    AddVessel(p_new_vessel);

    return p_new_vessel;

}

#ifdef CHASTE_VTK
template <unsigned DIM>
vtkSmartPointer<vtkPolyData> CaVascularNetwork<DIM>::GetVtk()
{
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
    unsigned numberOfNodes = 0;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it2;
    for(it2 = mVessels.begin(); it2 < mVessels.end(); it2++)
    {
        numberOfNodes += (*it2)->GetNumberOfNodes();
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
    unsigned node_index = 0;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it < mVessels.end(); it++)
    {
        vtkSmartPointer<vtkLine> pLine = vtkSmartPointer<vtkLine>::New();
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it)->GetSegments();

        for(unsigned i = 0; i < segments.size(); i++)
        {
            // Add the point for each node
            ChastePoint<DIM> location = segments[i]->GetNode(0)->GetLocation();
            vtkIdType pointId;
            if(DIM == 2)
            {
                pointId = pPoints->InsertNextPoint(location[0], location[1], 0.0);
            }
            else
            {
                pointId = pPoints->InsertNextPoint(location[0], location[1], location[2]);
            }

            std::map<std::string, double> vtk_node_data = segments[i]->GetNode(0)->GetVtkData();
            std::map<std::string, boost::any> generic_node_data = segments[i]->GetNode(0)->rGetDataContainer().GetMap();
            // Add the node data
            for(unsigned idx=0; idx < pNodeInfoVector.size(); idx++)
            {
                // Get the key
                std::string key = pNodeInfoVector[idx]->GetName();

                // If it is in the vtk data use it
                if(vtk_node_data.count(key) == 1)
                {
                    pNodeInfoVector[idx]->SetValue(node_index, vtk_node_data[key]);
                }
                // Otherwise check the generic data
                else if(generic_node_data.count(key) == 1)
                {
                    if(generic_node_data[key].type() == typeid(double))
                    {
                        double cast_value = boost::any_cast<double>(generic_node_data[key]);
                        pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                    }
                    else if(generic_node_data[key].type() == typeid(unsigned))
                    {
                        double cast_value = double(boost::any_cast<unsigned>(generic_node_data[key]));
                        pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                    }
                    else if(generic_node_data[key].type() == typeid(bool))
                    {
                        double cast_value = double(boost::any_cast<bool>(generic_node_data[key]));
                        pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                    }
                }
            }
            node_index++;
            pLine->GetPointIds()->InsertId(i, pointId);

            // Do an extra insert for the last node in the segment
            if (i == segments.size() - 1)
            {
                // Add the point for each node
                ChastePoint<DIM> location = segments[i]->GetNode(1)->GetLocation();
                vtkIdType pointId2;
                if(DIM == 2)
                {
                    pointId2 = pPoints->InsertNextPoint(location[0], location[1], 0.0);
                }
                else
                {
                    pointId2 = pPoints->InsertNextPoint(location[0], location[1], location[2]);
                }

                std::map<std::string, double> vtk_node_data = segments[i]->GetNode(1)->GetVtkData();
                std::map<std::string, boost::any> generic_node_data = segments[i]->GetNode(1)->rGetDataContainer().GetMap();
                // Add the node data
                for(unsigned idx=0; idx < pNodeInfoVector.size(); idx++)
                {
                    // Get the key
                    std::string key = pNodeInfoVector[idx]->GetName();

                    // If it is in the vtk data use it
                    if(vtk_node_data.count(key) == 1)
                    {
                        pNodeInfoVector[idx]->SetValue(node_index, vtk_node_data[key]);
                    }
                    // Otherwise check the generic data
                    else if(generic_node_data.count(key) == 1)
                    {
                        if(generic_node_data[key].type() == typeid(double))
                        {
                            double cast_value = boost::any_cast<double>(generic_node_data[key]);
                            pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                        }
                        else if(generic_node_data[key].type() == typeid(unsigned))
                        {
                            double cast_value = double(boost::any_cast<unsigned>(generic_node_data[key]));
                            pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                        }
                        else if(generic_node_data[key].type() == typeid(bool))
                        {
                            double cast_value = double(boost::any_cast<bool>(generic_node_data[key]));
                            pNodeInfoVector[idx]->SetValue(node_index, cast_value);
                        }
                    }
                }
                node_index++;
                pLine->GetPointIds()->InsertId(i + 1, pointId2);
            }
        }
        pLines->InsertNextCell(pLine);

        // Add the vessel data
        std::map<std::string, double> vtk_vessel_data = (*it)->GetVtkData();
        std::map<std::string, boost::any> generic_vessel_data = (*it)->rGetDataContainer().GetMap();
        // Add the node data
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
void CaVascularNetwork<DIM>::Write(const std::string& filename)
{
    vtkSmartPointer<vtkPolyData> p_polydata = GetVtk();
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(p_polydata);
    writer->Write();
}
#endif // CHASTE_VTK

template <unsigned DIM>
void CaVascularNetwork<DIM>::WriteConnectivity(const std::string& output_filename)
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

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iterator;
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


template<unsigned DIM>
void CaVascularNetwork<DIM>::UpdateNodes()
{
    mNodes = std::vector<boost::shared_ptr<VascularNode<DIM> > >();
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        std::vector<boost::shared_ptr<VascularNode<DIM> > > vessel_nodes = (*it)->GetNodes();
        std::copy(vessel_nodes.begin(), vessel_nodes.end(), std::inserter(nodes, nodes.begin()));
    }
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(mNodes));
    mNodesUpToDate = true;
}

template<unsigned DIM>
void CaVascularNetwork<DIM>::UpdateSegments()
{
    mSegments = std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >();
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > vessel_segments = (*it)->GetSegments();
        std::copy(vessel_segments.begin(), vessel_segments.end(), std::back_inserter(mSegments));
    }
    mSegmentsUpToDate = true;
}

template<unsigned DIM>
void CaVascularNetwork<DIM>::UpdateVesselNodes()
{
    mVesselNodes = std::vector<boost::shared_ptr<VascularNode<DIM> > >();
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        boost::shared_ptr<VascularNode<DIM> > vessel_node1 = (*it)->GetStartNode();
        nodes.insert(vessel_node1);
        boost::shared_ptr<VascularNode<DIM> > vessel_node2 = (*it)->GetEndNode();
        nodes.insert(vessel_node2);
    }
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(mVesselNodes));
    mVesselNodesUpToDate = true;
}

// Explicit instantiation
template class CaVascularNetwork<2>;
template class CaVascularNetwork<3>;

