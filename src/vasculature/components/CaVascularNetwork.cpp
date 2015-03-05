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

#include "CaVascularNetwork.hpp"

#include <boost/graph/graphviz.hpp>
#include <algorithm>
#include "OutputFileHandler.hpp"

template <unsigned DIM>
CaVascularNetwork<DIM>::CaVascularNetwork()
: mVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > >()),
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
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels)
{
    mVessels.insert(mVessels.end(), vessels.begin(), vessels.end());
}

template <unsigned DIM>
std::vector<std::pair<double, double> > CaVascularNetwork<DIM>::GetExtents()
{
    double x_max = -DBL_MAX;
    double y_max = -DBL_MAX;
    double z_max = -DBL_MAX;

    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it;

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
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::GetNearestNode(ChastePoint<DIM>& rLocation)
{
    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    boost::shared_ptr<VascularNode<DIM> > nearest_node;
    double min_distance = 1.e12;

    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {

        double node_distance = (*node_iter)->GetDistance(rLocation);
        if (node_distance < min_distance)
        {
            min_distance = node_distance;
            nearest_node = (*node_iter) ;
        }
    }
    return nearest_node;
}

template <unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVascularNetwork<DIM>::GetNearestSegment(const ChastePoint<DIM>& rLocation)
{
    boost::shared_ptr<CaVesselSegment<DIM> > nearest_segment;
    double min_distance = 1.e12;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iter;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter;
    for(vessel_iter = mVessels.begin(); vessel_iter != mVessels.end(); vessel_iter++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*vessel_iter)->GetSegments();
        for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
        {
            double segment_distance = (*segment_iter)->GetDistance(rLocation);
            if (segment_distance < min_distance)
            {
                min_distance = segment_distance;
                nearest_segment = (*segment_iter) ;
            }
        }
    }
    return nearest_segment;
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetNearestVessel(const ChastePoint<DIM>& rLocation)
{
    return GetNearestSegment(rLocation)->GetVessel();
}

template <unsigned DIM>
std::set<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::GetNodes()
{
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator jt;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it)->GetSegments();
        for(jt = segments.begin(); jt != segments.end(); jt++)
        {
            nodes.insert((*jt)->GetNode(0));
            nodes.insert((*jt)->GetNode(1));
        }
    }
    return nodes;
}

template <unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVascularNetwork<DIM>::GetNode(unsigned index)
{

    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes = GetNodes();

    std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it = nodes.begin();
    std::advance(it, index);

    return (*it);

}

template <unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::GetVectorOfNodes()
{
	std::vector<boost::shared_ptr<VascularNode<DIM> > > vec_nodes;
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes = GetNodes();
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(vec_nodes));

    return vec_nodes;
}


template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfNodes()
{
    return GetNodes().size();
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
std::set<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::GetVesselEndNodes()
{
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator jt;

    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        nodes.insert((*it)->GetEndNode());
        nodes.insert((*it)->GetStartNode());
    }
    return nodes;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::GetVessels()
{
    return mVessels;
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
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;

    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;

    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        for(unsigned i = 0; i != (*it)->GetNumberOfSegments(); i++)
        {
            segments.push_back((*it)->GetSegment(i));
        }
    }

    return segments;
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
    std::set<boost::shared_ptr<VascularNode<DIM> > >  vessel_nodes = GetVesselEndNodes();
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
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
std::vector<std::vector<bool> > CaVascularNetwork<DIM>::IsConnected(std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes,
                                                                    std::vector<boost::shared_ptr<VascularNode<DIM> > > queryNodes)
                                                                    {

    // Assign the vessel nodes unique IDs
    std::set<boost::shared_ptr<VascularNode<DIM> > >  vessel_nodes = GetVesselEndNodes();
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
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

    std::vector<std::vector<bool> > connected;

    for(unsigned i=0; i<sourceNodes.size(); i++)
    {
        boost::shared_ptr<VascularNode<DIM> > pSourceNode = sourceNodes[i];
        boost::shared_ptr<CaVessel<DIM> > p_source_vessel = pSourceNode->GetVesselSegment(0)->GetVessel();
        boost::shared_ptr<VascularNode<DIM> > pEquivalentSourceNode = p_source_vessel->GetStartNode();

        std::vector<bool> connected_to_this_source;

        for (unsigned j=0; j<queryNodes.size(); j++)
        {
            boost::shared_ptr<VascularNode<DIM> > pQueryNode = queryNodes[j];

            if (pSourceNode == pQueryNode)
            {
                connected_to_this_source.push_back(true);
                continue;
            }

            boost::shared_ptr<CaVessel<DIM> > p_query_vessel = pQueryNode->GetVesselSegment(0)->GetVessel();
            if (p_source_vessel == p_query_vessel || p_source_vessel->IsConnectedTo(p_query_vessel))
            {
                connected_to_this_source.push_back(true);
                continue;
            }

            boost::shared_ptr<VascularNode<DIM> > pEquivalentQueryNode = p_query_vessel->GetStartNode();


            // a vector to hold the discover time property for each vertex
            std::vector<Size> dtime(num_vertices(G));

            Size time = 0;
            bfs_time_visitor<Size*>vis(&dtime[0], time);

            breadth_first_search(G,vertex(pEquivalentSourceNode->GetId(),G), boost::visitor(vis));

            connected_to_this_source.push_back(dtime[pEquivalentQueryNode->GetId()] > 0);
        }

        connected.push_back(connected_to_this_source);
    }

    return connected;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::MergeCoincidentNodes()
{
    // Loop through the nodes, if they are coincident but not identical replace
    // one of them. Nodes in the inner iteration loop are the ones replaced.
    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();

    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it;
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it2;
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
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetNodeData(VasculatureData data)
{
    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = GetNodes();
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetVesselData(VasculatureData data)
{
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it != mVessels.end(); it++)
    {
        (*it)->SetDataContainer(data);
    }
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::Translate(const std::vector<double>& rTranslationVector, bool copy)
{
    //	Get a vector of points corresponding to nodes, move the points
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;
    std::vector<ChastePoint<DIM> > points;

    if(!copy)
    {
        nodes = GetNodes();
    }
    else
    {
        nodes = DeepCopyVessels();
    }

    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {
        points.push_back((*node_iter)->GetLocation());
    }

    GeometryTransform<DIM> transform;
    std::vector<ChastePoint<DIM> > new_points = transform.Translate(points, rTranslationVector);

    unsigned count = 0u;
    for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
    {
        (*node_iter)->SetLocation(new_points[count]);
        count++;
    }

    // Merge coincident nodes
    if(copy)
    {
        MergeCoincidentNodes();
    }
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::Write(const std::string& filename, bool geometry_only)
{
    vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pPoints= vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> pLines = vtkSmartPointer<vtkCellArray>::New();

    std::vector<vtkSmartPointer<vtkDoubleArray> > pVesselInfoVector;
    std::vector<vtkSmartPointer<vtkDoubleArray> > pNodeInfoVector;

    // Set up vessel info arrays
    std::map<std::string, boost::any >::iterator map_iterator;
    std::map<std::string, boost::any > data_map = mVessels[0]->rGetDataContainer().GetMap();

    for(map_iterator = data_map.begin(); map_iterator != data_map.end(); map_iterator++)
    {
        vtkSmartPointer<vtkDoubleArray> pVesselInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pVesselInfo->SetNumberOfComponents(1); // all scalar data - has one entry
        pVesselInfo->SetNumberOfTuples(mVessels.size()); // number of tuples is number of vessels
        pVesselInfo->SetName((*map_iterator).first.c_str());

        // Only write information that can be cast to double
        if (map_iterator->second.type() == typeid(double) || map_iterator->second.type() == typeid(unsigned) ||
                map_iterator->second.type() == typeid(bool))
        {
            pVesselInfoVector.push_back(pVesselInfo);
        }
    }

    unsigned numberOfNodes = 0;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it2;
    for(it2 = mVessels.begin(); it2 < mVessels.end(); it2++)
    {
        numberOfNodes += (*it2)->GetNumberOfNodes();
    }

    // Set up node info arrays
    std::map<std::string, boost::any >::iterator map_iterator5;
    std::map<std::string, boost::any > data_map2 = mVessels[0]->GetStartNode()->rGetDataContainer().GetMap();
    for(map_iterator5 = data_map2.begin(); map_iterator5 != data_map2.end(); map_iterator5++)
    {
        vtkSmartPointer<vtkDoubleArray> pNodeInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pNodeInfo->SetNumberOfComponents(1); // all scalar data - has one entry
        pNodeInfo->SetNumberOfTuples(numberOfNodes); // number of tuples is number of vessels
        pNodeInfo->SetName((*map_iterator5).first.c_str());

        // Only write information that can be cast to double
        if (map_iterator5->second.type() == typeid(double) || map_iterator5->second.type() == typeid(unsigned) ||
                map_iterator5->second.type() == typeid(bool))
        {
            pNodeInfoVector.push_back(pNodeInfo);
        }
    }

    unsigned vessel_index=0;
    unsigned node_index = 0;
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
    for(it = mVessels.begin(); it < mVessels.end(); it++)
    {
        vtkSmartPointer<vtkLine> pLine = vtkSmartPointer<vtkLine>::New();
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it)->GetSegments();


        for(unsigned i = 0; i < segments.size(); i++)
        {
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
            std::map<std::string, boost::any >::iterator map_iterator3;
            unsigned key_index = 0;

            std::map<std::string, boost::any > data_map3 = segments[i]->GetNode(0)->rGetDataContainer().GetMap();
            for(map_iterator3 = data_map3.begin(); map_iterator3 != data_map3.end(); map_iterator3++)
            {
                if (map_iterator3->second.type() == typeid(double))
                {
                    double cast_value = boost::any_cast<double>(map_iterator3->second);
                    pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                    key_index++;
                }
                if (map_iterator3->second.type() == typeid(unsigned))
                {
                    double cast_value = (double)boost::any_cast<unsigned>(map_iterator3->second);
                    pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                    key_index++;
                }
                if (map_iterator3->second.type() == typeid(bool))
                {
                    double cast_value = (double)boost::any_cast<bool>(map_iterator3->second);
                    pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                    key_index++;
                }
            }
            node_index++;

            pLine->GetPointIds()->InsertId(i, pointId);
            if (i == segments.size() - 1)
            {
                vtkIdType pointId2;
                ChastePoint<DIM> location2 = segments[i]->GetNode(1)->GetLocation();
                if(DIM == 2)
                {
                    pointId2 = pPoints->InsertNextPoint(location2[0], location2[1], 0.0);
                }
                else
                {
                    pointId2 = pPoints->InsertNextPoint(location2[0], location2[1], location2[2]);
                }
                pLine->GetPointIds()->InsertId(i + 1, pointId2);
                std::map<std::string, boost::any >::iterator map_iterator4;

                unsigned key_index = 0;
                std::map<std::string, boost::any > data_map4 = segments[i]->GetNode(1)->rGetDataContainer().GetMap();
                for(map_iterator4 = data_map4.begin(); map_iterator4 != data_map4.end(); map_iterator4++)
                {
                    if (map_iterator4->second.type() == typeid(double))
                    {
                        double cast_value = boost::any_cast<double>(map_iterator4->second);
                        pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                        key_index++;
                    }
                    if (map_iterator4->second.type() == typeid(unsigned))
                    {
                        double cast_value = (double)boost::any_cast<unsigned>(map_iterator4->second);
                        pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                        key_index++;
                    }
                    if (map_iterator4->second.type() == typeid(bool))
                    {
                        double cast_value = (double)boost::any_cast<bool>(map_iterator4->second);
                        pNodeInfoVector[key_index]->SetValue(node_index, cast_value );
                        key_index++;
                    }
                }
                node_index++;
            }
        }
        pLines->InsertNextCell(pLine);

        std::map<std::string, boost::any >::iterator map_iterator2;
        std::map<std::string, boost::any > data_map5 = (*it)->rGetDataContainer().GetMap();

        unsigned key_index = 0;
        for(map_iterator2 = data_map5.begin(); map_iterator2 != data_map5.end(); map_iterator2++)
        {
            if (map_iterator2->second.type() == typeid(double))
            {
                double cast_value = boost::any_cast<double>(map_iterator2->second);
                pVesselInfoVector[key_index]->SetValue(vessel_index, cast_value );
                key_index++;
            }
            if (map_iterator2->second.type() == typeid(unsigned))
            {
                double cast_value = (double)boost::any_cast<unsigned>(map_iterator2->second);
                pVesselInfoVector[key_index]->SetValue(vessel_index, cast_value );
                key_index++;
            }
            if (map_iterator2->second.type() == typeid(bool))
            {
                double cast_value = (double)boost::any_cast<bool>(map_iterator2->second);
                pVesselInfoVector[key_index]->SetValue(vessel_index, cast_value );
                key_index++;
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

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(pPolyData);
    writer->Write();
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::VisualiseVesselConnectivity(std::string output_filename)
{

    // construct graph representation of vessel network
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

    Graph G;

    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iterator;

    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = GetVesselEndNodes();

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

template <unsigned DIM>
std::set<boost::shared_ptr<VascularNode<DIM> > > CaVascularNetwork<DIM>::DeepCopyVessels()
{
    std::set<boost::shared_ptr<VascularNode<DIM> > >  nodes;

    // Deep copy each of the vessels
    // loop through the vessels, segments and nodes and recreate new instances
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iter;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter;
    std::vector<boost::shared_ptr<CaVessel<DIM> > > new_vessels;

    for(vessel_iter = mVessels.begin(); vessel_iter != mVessels.end(); vessel_iter++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > new_segments;
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*vessel_iter)->GetSegments();

        for(segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
        {
            ChastePoint<DIM> node0_location = (*segment_iter)->GetNode(0)->GetLocation();
            ChastePoint<DIM> node1_location = (*segment_iter)->GetNode(1)->GetLocation();

            boost::shared_ptr<VascularNode<DIM> > pNode0(VascularNode<DIM>::Create(node0_location));
            boost::shared_ptr<VascularNode<DIM> > pNode1(VascularNode<DIM>::Create(node1_location));
            nodes.insert(pNode0);
            nodes.insert(pNode1);

            boost::shared_ptr<CaVesselSegment<DIM> > pSegment(CaVesselSegment<DIM>::Create(pNode0, pNode1));
            new_segments.push_back(pSegment);
        }
        boost::shared_ptr<CaVessel<DIM> > pVessel(CaVessel<DIM>::Create(new_segments));
        new_vessels.push_back(pVessel);
    }
    AddVessels(new_vessels);

    return nodes;
}

// Explicit instantiation
template class CaVascularNetwork<2>;
template class CaVascularNetwork<3>;

