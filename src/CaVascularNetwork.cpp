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

template <unsigned DIM>
CaVascularNetwork<DIM>::CaVascularNetwork()
	: mVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > >()),
	  //mNodes(std::vector<boost::shared_ptr<CaVascularNetworkNode<DIM> > >()),
	  mpDataContainer(boost::shared_ptr<VascularNetworkData>(new VascularNetworkData()))
{
}

template <unsigned DIM>
CaVascularNetwork<DIM>::~CaVascularNetwork()
{
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessels(boost::shared_ptr<CaVessel<DIM> > vessel)
{
	mVessels.push_back(vessel);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels)
{
	mVessels.insert(mVessels.end(), vessels.begin(), vessels.end());
}

template <unsigned DIM>
std::set<boost::shared_ptr<CaVascularNetworkNode<DIM> > > CaVascularNetwork<DIM>::GetNodes()
{
	std::set<boost::shared_ptr<CaVascularNetworkNode<DIM> > >  nodes;

	typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;
	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator jt;
	for(it = mVessels.begin(); it != mVessels.end(); it++)
	{
		std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it)->GetSegments();
		for(jt = segments.begin(); jt != segments.end(); jt++)
		{
			nodes.insert((*jt)->GetNodes(0));
			nodes.insert((*jt)->GetNodes(1));
		}
	}
	return nodes;
}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::GetVessels()
{
	return mVessels;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::MergeCoincidentNodes()
{
	// Loop through the nodes, if they are coincident but not identical replace
	// one of them. Nodes in the inner iteration loop are the ones replaced.
	std::set<boost::shared_ptr<CaVascularNetworkNode<DIM> > > nodes = GetNodes();

	typename std::set<boost::shared_ptr<CaVascularNetworkNode<DIM> > >::iterator it;
	typename std::set<boost::shared_ptr<CaVascularNetworkNode<DIM> > >::iterator it2;
	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it3;

	for(it = nodes.begin(); it != nodes.end(); it++)
	{
		for(it2 = nodes.begin(); it2 != nodes.end(); it2++)
		{
			// If the nodes are not identical
			if ((*it) != (*it2))
			{
				// If the node locations are the same - according to the ChastePoint definition
				if((*it)->GetLocation().IsSamePoint((*it2)->GetLocation()))
				{
					// Replace the node corresponding to it2 with the one corresponding to it
					// in all segments.
					std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it2)->GetVesselSegments();
					for(it3 = segments.begin(); it3 != segments.end(); it3++)
					{
						if ((*it3)->GetNodes(0) == (*it2))
						{
							(*it3)->ReplaceNode(0, (*it));
						}
						else if(((*it3)->GetNodes(1) == (*it2)))
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
void CaVascularNetwork<DIM>::WriteToFile(std::string filename, bool geometry_only)
{
	vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pPoints= vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> pLines = vtkSmartPointer<vtkCellArray>::New();
	//vtkSmartPointer<vtkFloatArray> pInfo = vtkSmartPointer<vtkFloatArray>::New();

	//pInfo->SetNumberOfComponents(1);
	//pInfo->SetName("Something");

	typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator it;

    for(it = mVessels.begin(); it < mVessels.end(); it++)
    {
    	vtkSmartPointer<vtkLine> pLine = vtkSmartPointer<vtkLine>::New();
    	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*it)->GetSegments();

    	for(unsigned i = 0; i < segments.size(); i++)
    	{
    		ChastePoint<DIM> location = segments[i]->GetNodes(0)->GetLocation();
    		vtkIdType pointId;
    		if(DIM == 2)
    		{
    			pointId = pPoints->InsertNextPoint(location[0], location[1], 0.0);
    		}
    		else
    		{
    			pointId = pPoints->InsertNextPoint(location[0], location[1], location[2]);
    		}
    		pLine->GetPointIds()->InsertId(i, pointId);

    		if (i == segments.size() - 1)
    		{
    			vtkIdType pointId2;
    			ChastePoint<DIM> location2 = segments[i]->GetNodes(1)->GetLocation();
        		if(DIM == 2)
        		{
        			pointId2 = pPoints->InsertNextPoint(location2[0], location2[1], 0.0);
        		}
        		else
        		{
        			pointId2 = pPoints->InsertNextPoint(location2[0], location2[1], location2[2]);
        		}
        		pLine->GetPointIds()->InsertId(i + 1, pointId2);
    		}
    	}
    	pLines->InsertNextCell(pLine);
    	//pInfo->InsertNextTupleValue(1.0);
    }
    pPolyData->SetPoints(pPoints);
    pPolyData->SetLines(pLines);
    //pPolyData->GetCellData().SetScalars(pInfo);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInput(pPolyData);
	writer->Write();
}


//template <unsigned DIM>
//boost::shared_ptr<CaVascularNetwork<DIM> > CaVascularNetwork<DIM>::shared()
//{
//    return this->shared_from_this();
//}
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::GetVesselID(boost::shared_ptr<CaVessel<DIM> > vessel)
//{
//    bool vessel_found = false;
//    unsigned i = 0;
//
//    for (i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (vessel == mVesselArray[i])
//        {
//        	vessel_found = true;
//            break;
//        }
//    }
//
//    assert(vessel_found);
//    return i;
//}
//
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::GetNodeID(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
//{
//    bool node_found = false;
//    unsigned i = 0;
//
//    for (i = 0; i < GetNumberOfNodesInNetwork(); i++)
//    {
//        if (node == mNodeArray[i])
//        {
//        	node_found = true;
//            break;
//        }
//    }
//
//    assert(node_found);
//    return i;
//}
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::GetNumberOfVesselsInNetwork()
//{
//    return mVesselArray.size();
//}
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::GetNumberOfNodesInNetwork()
//{
//    return mNodeArray.size();
//}
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::GetNumberOfVesselsAtLocation(ChastePoint<DIM> coord)
//{
//
//	unsigned numberOfVesselsAtLocation = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
//        {
//            if (GetVessel(i)->GetSegmentCoordinate(j).IsSamePoint(coord))
//            {
//                numberOfVesselsAtLocation++;
//            }
//        }
//    }
//
//    return numberOfVesselsAtLocation;
//}
//
//
//template <unsigned DIM>
//bool CaVascularNetwork<DIM>::NodePresentAtLocation(ChastePoint<DIM> location)
//{
//    bool nodePresentAtLocation = false;
//
//    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
//    {
//        if (mNodeArray[i]->GetLocation().IsSamePoint(location))
//        {
//            nodePresentAtLocation = true;
//            break;
//        }
//    }
//    return nodePresentAtLocation;
//}
//
//template <unsigned DIM>
//unsigned CaVascularNetwork<DIM>::NumberOfNodesPresentAtLocation(ChastePoint<DIM> location)
//{
//	unsigned numberOfNodesPresentAtLocation = 0;
//
//    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
//    {
//        if (mNodeArray[i]->GetLocation().IsSamePoint(location))
//        {
//            numberOfNodesPresentAtLocation++;
//        }
//    }
//
//    return numberOfNodesPresentAtLocation;
//}
//
//template <unsigned DIM>
//bool CaVascularNetwork<DIM>::VesselIsInNetwork(boost::shared_ptr<CaVessel<DIM> > vessel)
//{
//    bool vesselIsInNetwork = false;
//    unsigned i = 0;
//
//    for (i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (vessel == mVesselArray[i])
//        {
//            vesselIsInNetwork = true;
//            break;
//        }
//    }
//    return vesselIsInNetwork;
//}
//
//template <unsigned DIM>
//bool CaVascularNetwork<DIM>::NodeIsInNetwork(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
//{
//    bool nodeIsInNetwork = false;
//
//    unsigned i = 0;
//
//    for (i = 0; i < GetNumberOfNodesInNetwork(); i++)
//    {
//        if (node == mNodeArray[i])
//        {
//            nodeIsInNetwork = true;
//            break;
//        }
//    }
//
//    return nodeIsInNetwork;
//}
//
///*
// * Helper class for "connected" methods
// */
//template<typename TimeMap> class bfs_time_visitor : public boost::default_bfs_visitor
//{
//    typedef typename boost::property_traits<TimeMap>::value_type T;
//
//public:
//
//    TimeMap m_timemap;
//    T& m_time;
//
//    bfs_time_visitor(TimeMap tmap, T& t)
//    	:m_timemap(tmap),
//    	 m_time(t)
//    {
//    }
//
//    template<typename Vertex, typename Graph>
//    void discover_vertex(Vertex u, const Graph& g) const
//    {
//        put(m_timemap, u, m_time++);
//    }
//};
//
/////\ todo this could be made more general by passing in vectors of source nodes and target nodes and returning true for any targets connected
//// to a source. Would avoid graph reconstruction for each query and results in only one overall method.
//template <unsigned DIM>
//bool CaVascularNetwork<DIM>::Connected(boost::shared_ptr<CaVascularNetworkNode<DIM> > node1, boost::shared_ptr<CaVascularNetworkNode<DIM> > node2)
//{
//
//    if (node1 == node2)
//    {
//        return true;
//    }
//
//    // construct graph representation of vessel network
//    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
//
//    Graph G;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        add_edge(GetNodeID(GetVessel(i)->GetNode1()), GetNodeID(GetVessel(i)->GetNode2()), G);
//    }
//
//    // typedefs
//    typedef boost::graph_traits<Graph>::vertices_size_type Size;
//
//    // a vector to hold the discover time property for each vertex
//    std::vector<Size> dtime(num_vertices(G));
//
//    Size time = 0;
//    bfs_time_visitor<Size*>vis(&dtime[0], time);
//
//    // use breadth first search to establish discovery time of all nodes from node1
//    // this assigns a discovery time to dTime for each node (index relates to nodeID)
//    // dTime is zero for node1 and all other nodes that are not connected to node1
//    // dTime is nonzero for all nodes that are connected to node1 (except node 1 itself)
//    breadth_first_search(G,vertex(GetNodeID(node1),G), boost::visitor(vis));
//
//    return (dtime[GetNodeID(node2)] > 0);
//}
//
//template <unsigned DIM>
//bool CaVascularNetwork<DIM>::ConnectedToInputNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
//{
//
//    if (node->GetBooleanData("IsInputNode"))
//    {
//        return true;
//    }
//
//    // construct graph representation of vessel network
//    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;
//
//    Graph G;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        add_edge(GetNodeID(GetVessel(i)->GetNode1()),GetNodeID(GetVessel(i)->GetNode2()), G);
//    }
//
//    // typedefs
//    typedef boost::graph_traits<Graph>::vertices_size_type Size;
//
//    int connectedToInputNode = 0;
//
//    for (unsigned j = 0; j < GetNumberOfNodesInNetwork(); j++)
//    {
//
//        if (GetNode(j)->GetBooleanData("IsInputNode"))
//        {
//            // a vector to hold the discover time property for each vertex
//            std::vector<Size> dtime(num_vertices(G));
//
//            Size time = 0;
//            bfs_time_visitor<Size*>vis(&dtime[0], time);
//
//            // use breadth first search to establish discovery time of all nodes from node1
//            // this assigns a discovery time to dTime for each node (index relates to nodeID)
//            // dTime is zero for node1 and all other nodes that are not connected to node1
//            // dTime is nonzero for all nodes that are connected to node1 (except node 1 itself)
//            breadth_first_search(G,vertex(j,G), boost::visitor(vis));
//
//            if (dtime[GetNodeID(node)] > 0)
//            {
//                connectedToInputNode++;
//            }
//        }
//    }
//    return (connectedToInputNode > 0);
//}
//

// Explicit instantiation
template class CaVascularNetwork<2>;
template class CaVascularNetwork<3>;

