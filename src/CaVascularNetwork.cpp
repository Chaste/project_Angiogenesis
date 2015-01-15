/*
 * CaVascularNetwork.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#include "CaVascularNetwork.hpp"

//
//  CaVascularNetwork.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 16/10/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include "CaVascularNetwork.hpp"
#include <math.h>
#include <float.h>

#include <boost/config.hpp>
#include <algorithm>
#include <utility>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <boost/pending/indirect_cmp.hpp>



template <unsigned SPATIAL_DIM>
CaVascularNetwork<SPATIAL_DIM>::CaVascularNetwork() : mVesselArray(), mNodeArray(), mArterialHaematocritLevel(0.45), mArterialInputPressure(27*(1.01*pow(10.0,5)/760)), mVenousOutputPressure(15*(1.01*pow(10.0,5)/760))
{

}

template <unsigned SPATIAL_DIM>
CaVascularNetwork<SPATIAL_DIM>::~CaVascularNetwork()
{

}

template <unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > CaVascularNetwork<SPATIAL_DIM>::shared()
{
    return this->shared_from_this();
}

template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::GetVesselID(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{

    bool vesselFound = false;
    unsigned i = 0;

    for (i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (vessel == mVesselArray[i])
        {
            vesselFound = true;
            break;
        }
    }

    assert(vesselFound);

    return i;
}


template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::GetNodeID(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{

    bool nodeFound = false;
    unsigned i = 0;

    for (i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (node == mNodeArray[i])
        {
            nodeFound = true;
            break;
        }
    }

    assert(nodeFound);

    return i;

}

template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::GetNumberOfVesselsInNetwork()
{
    return mVesselArray.size();
}

template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::GetNumberOfNodesInNetwork()
{
    return mNodeArray.size();
}


template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::GetNumberOfVesselsAtLocation(ChastePoint<SPATIAL_DIM> coord)
{

	unsigned numberOfVesselsAtLocation = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
        {
            if (GetVessel(i)->GetSegmentCoordinate(j).IsSamePoint(coord))
            {
                numberOfVesselsAtLocation++;
            }
        }
    }

    return numberOfVesselsAtLocation;
}

template <unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > CaVascularNetwork<SPATIAL_DIM>::GetVessel(int vessel_id)
{
    return mVesselArray[vessel_id];
}

template <unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > CaVascularNetwork<SPATIAL_DIM>::GetVessel(ChastePoint<SPATIAL_DIM> coord, int positionInContainer)
{

    int numberOfVesselsAtLocation = 0;
    int vesselID = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
        {
            if (GetVessel(i)->GetSegmentCoordinate(j).IsSamePoint(coord))
            {
                numberOfVesselsAtLocation++;
            }

            if (numberOfVesselsAtLocation > positionInContainer)
            {
                vesselID = i;
            }
        }
    }

    if (numberOfVesselsAtLocation <= positionInContainer)
    {
        // todo This should strictly be an exception.
        assert(numberOfVesselsAtLocation > positionInContainer);
    }

    return mVesselArray[vesselID];
}

template <unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > CaVascularNetwork<SPATIAL_DIM>::GetNode(int node_id)
{
    return mNodeArray[node_id];
}

template <unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > CaVascularNetwork<SPATIAL_DIM>::GetNode(ChastePoint<SPATIAL_DIM> location)
{

    assert(NodePresentAtLocation(location));
    assert(NumberOfNodesPresentAtLocation(location) == 1);

    int nodeID;
    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (mNodeArray[i]->GetLocation().IsSamePoint(location))
        {
            nodeID = i;
            break;
        }
    }

    return mNodeArray[nodeID];

}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetArterialHaematocritLevel()
{
    return mArterialHaematocritLevel;
}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetArterialInputPressure()
{
    return mArterialInputPressure;
}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetVenousOutputPressure()
{
    return mVenousOutputPressure;
}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetMeanVesselLengthOfNeovasculature()
{
    double totalVesselLength = 0;
    int numberofNewVessels = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (GetVessel(i)->IsPartOfNeovasculature())
        {
            totalVesselLength += GetVessel(i)->GetLength();
            numberofNewVessels++;
        }

    }

    return (totalVesselLength/(double)numberofNewVessels);
}

template <unsigned SPATIAL_DIM>
int CaVascularNetwork<SPATIAL_DIM>::GetNumberOfVesselsByLength(double lowerBoundLength, double upperBoundLength)
{
    int numberOfVesselsInsideRange = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (GetVessel(i)->GetLength() > lowerBoundLength && GetVessel(i)->GetLength() <= upperBoundLength)
        {
            numberOfVesselsInsideRange++;
        }
    }

    return numberOfVesselsInsideRange;
}

template <unsigned SPATIAL_DIM>
int CaVascularNetwork<SPATIAL_DIM>::GetNumberOfVesselsByRadius(double lowerBoundRadius, double upperBoundRadius)
{
    int numberOfVesselsInsideRange = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (GetVessel(i)->GetRadius() > lowerBoundRadius && GetVessel(i)->GetRadius() <= upperBoundRadius)
        {
            numberOfVesselsInsideRange++;
        }
    }

    return numberOfVesselsInsideRange;
}

template <unsigned SPATIAL_DIM>
int CaVascularNetwork<SPATIAL_DIM>::GetNumberOfVesselsByTortuosity(double lowerBoundTortuosity, double upperBoundTortuosity)
{
    int numberOfVesselsInsideRange = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (GetVessel(i)->GetTortuosity() > lowerBoundTortuosity && GetVessel(i)->GetTortuosity() <= upperBoundTortuosity)
        {
            numberOfVesselsInsideRange++;
        }
    }

    return numberOfVesselsInsideRange;
}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetMeanVesselRadiusOfNeovasculature()
{
    double totalVesselRadius = 0;
    int numberofNewVessels = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (GetVessel(i)->IsPartOfNeovasculature())
        {
            totalVesselRadius += GetVessel(i)->GetRadius();
            numberofNewVessels++;
        }

    }

    return (totalVesselRadius/(double)numberofNewVessels);
}

template <unsigned SPATIAL_DIM>
double CaVascularNetwork<SPATIAL_DIM>::GetMeanVesselTortuosityOfNeovasculature()
{
    double totalVesselTortuosity = 0;
    int numberofNewVessels = 0;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (std::isinf(GetVessel(i)->GetTortuosity()))
        {
            // we ignore vessels with infinite tortuosity, i.e. self-loops
        }
        else
        {
            if (GetVessel(i)->IsPartOfNeovasculature())
            {
                totalVesselTortuosity += GetVessel(i)->GetTortuosity();
                numberofNewVessels++;
            }
        }
    }

    return (totalVesselTortuosity/(double)numberofNewVessels);
}

template <unsigned SPATIAL_DIM>
std::vector<boost::shared_ptr<CaVessel<SPATIAL_DIM> > > CaVascularNetwork<SPATIAL_DIM>::GetVessels()
{
    return mVesselArray;
}

template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SetArterialHaematocritLevel(double value)
{
    mArterialHaematocritLevel = value;
}

template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SetArterialInputPressure(double value)
{
    mArterialInputPressure = value;
}

template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SetVenousOutputPressure(double value)
{
    mVenousOutputPressure = value;
}

template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::NodePresentAtLocation(ChastePoint<SPATIAL_DIM> location)
{
    bool nodePresentAtLocation = false;

    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (mNodeArray[i]->GetLocation().IsSamePoint(location))
        {
            nodePresentAtLocation = true;
            break;
        }
    }

    return nodePresentAtLocation;
}

template <unsigned SPATIAL_DIM>
unsigned CaVascularNetwork<SPATIAL_DIM>::NumberOfNodesPresentAtLocation(ChastePoint<SPATIAL_DIM> location)
{
	unsigned numberOfNodesPresentAtLocation = 0;

    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (mNodeArray[i]->GetLocation().IsSamePoint(location))
        {
            numberOfNodesPresentAtLocation++;
        }
    }

    return numberOfNodesPresentAtLocation;
}

template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::AddVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{

    // todo Checking that there is enough space for the vessel at a particular location should be handled in a new CABasedCellPopulation

    // Check that the first and last segment coordinates are the same as the coordinates of the vessel nodes.

    assert(vessel->GetNode1()->GetLocation().IsSamePoint(vessel->GetSegmentCoordinate(0)) || vessel->GetNode2()->GetLocation().IsSamePoint(vessel->GetSegmentCoordinate(0)));
    assert(vessel->GetNode1()->GetLocation().IsSamePoint(vessel->GetSegmentCoordinate(vessel->GetNumberOfSegments() - 1)) || vessel->GetNode2()->GetLocation().IsSamePoint(vessel->GetSegmentCoordinate(vessel->GetNumberOfSegments() - 1)));

    assert(vessel->GetNode1()->GetNumberOfAdjoiningVessels() == 1);
    assert(vessel->GetNode2()->GetNumberOfAdjoiningVessels() == 1);

    // add vessel to network

    mVesselArray.push_back(vessel);

    bool node1OfVesselAlreadyPresentInNetwork = false;

    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (GetNode(i)->GetLocation().IsSamePoint(vessel->GetNode1()->GetLocation())) // assume nodes are the same if they have the same location
        {
            vessel->SetNode1(mNodeArray[i]);
            GetNode(i)->AddAdjoiningVessel(vessel);
            // we do not want to change the number of vessels which are attached to input or output
            // nodes of the network for the moment - this is a limitation of the models
            // some of the mathematics (haematocrit calculations) require that only one vessel exits/enters
            // an input/output node of the network currently
            // ----------------------------------------------------------------------------------
            assert(GetNode(i)->IsInputNode() == false);
            assert(GetNode(i)->IsOutputNode() == false);
            // ----------------------------------------------------------------------------------
            node1OfVesselAlreadyPresentInNetwork = true;
            break;
        }
    }

    if (node1OfVesselAlreadyPresentInNetwork == false)
    {
        mNodeArray.push_back(vessel->GetNode1());
    }

    bool node2OfVesselAlreadyPresentInNetwork = false;

    for (unsigned i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (GetNode(i)->GetLocation().IsSamePoint(vessel->GetNode2()->GetLocation())) // assume nodes are the same if they have the same location
        {
            vessel->SetNode2(mNodeArray[i]);
            GetNode(i)->AddAdjoiningVessel(vessel);
            // we do not want to change the number of vessels which are attached to input or output
            // nodes of the network for the moment - this is a limitation of the models
            // some of the mathematics (haematocrit calculations) require that only one vessel exits/enters
            // an input/output node of the network currently
            // ----------------------------------------------------------------------------------
            assert(GetNode(i)->IsInputNode() == false);
            assert(GetNode(i)->IsOutputNode() == false);
            // ----------------------------------------------------------------------------------
            node2OfVesselAlreadyPresentInNetwork = true;
            break;
        }
    }

    if (node2OfVesselAlreadyPresentInNetwork == false)
    {
        mNodeArray.push_back(vessel->GetNode2());
    }

}

template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::VesselIsInNetwork(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel)
{
    bool vesselIsInNetwork = false;

    unsigned i = 0;

    for (i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (vessel == mVesselArray[i])
        {
            vesselIsInNetwork = true;
            break;
        }
    }

    return vesselIsInNetwork;
}


template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::NodeIsInNetwork(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{
    bool nodeIsInNetwork = false;

    unsigned i = 0;

    for (i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (node == mNodeArray[i])
        {
            nodeIsInNetwork = true;
            break;
        }
    }

    return nodeIsInNetwork;
}

/*
 * Helper class for "connected" methods
 */
template < typename TimeMap > class bfs_time_visitor:public boost::default_bfs_visitor {
    typedef typename boost::property_traits < TimeMap >::value_type T;
public:
    bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
    {
        put(m_timemap, u, m_time++);
    }
    TimeMap m_timemap;
    T & m_time;
};

template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::Connected(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node1, boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node2)
{

    if (node1 == node2)
    {
        return true;
    }

    // construct graph representation of vessel network

    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(GetNodeID(GetVessel(i)->GetNode1()),GetNodeID(GetVessel(i)->GetNode2()),G);
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
    breadth_first_search(G,vertex(GetNodeID(node1),G), boost::visitor(vis));

    return (dtime[GetNodeID(node2)] > 0);

}

template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::ConnectedToInputNode(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{

    if (node->IsInputNode())
    {
        return true;
    }

    // construct graph representation of vessel network

    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(GetNodeID(GetVessel(i)->GetNode1()),GetNodeID(GetVessel(i)->GetNode2()),G);
    }

    // typedefs

    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    int connectedToInputNode = 0;

    for (unsigned j = 0; j < GetNumberOfNodesInNetwork(); j++)
    {

        if (GetNode(j)->IsInputNode())
        {
            // a vector to hold the discover time property for each vertex
            std::vector<Size> dtime(num_vertices(G));

            Size time = 0;
            bfs_time_visitor<Size*>vis(&dtime[0], time);

            // use breadth first search to establish discovery time of all nodes from node1
            // this assigns a discovery time to dTime for each node (index relates to nodeID)
            // dTime is zero for node1 and all other nodes that are not connected to node1
            // dTime is nonzero for all nodes that are connected to node1 (except node 1 itself)
            breadth_first_search(G,vertex(j,G), boost::visitor(vis));

            if (dtime[GetNodeID(node)] > 0)
            {
                connectedToInputNode++;
            }

        }

    }

    return (connectedToInputNode > 0);

}


template <unsigned SPATIAL_DIM>
bool CaVascularNetwork<SPATIAL_DIM>::ConnectedToOutputNode(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{

    if (node->IsOutputNode())
    {
        return true;
    }

    // construct graph representation of vessel network

    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(GetNodeID(GetVessel(i)->GetNode1()),GetNodeID(GetVessel(i)->GetNode2()),G);
    }

    // typedefs

    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    int connectedToOutputNode = 0;

    for (unsigned j = 0; j < GetNumberOfNodesInNetwork(); j++)
    {

        if (GetNode(j)->IsOutputNode())
        {
            // a vector to hold the discover time property for each vertex
            std::vector<Size> dtime(num_vertices(G));

            Size time = 0;
            bfs_time_visitor<Size*>vis(&dtime[0], time);

            // use breadth first search to establish discovery time of all nodes from node1
            // this assigns a discovery time to dTime for each node (index relates to nodeID)
            // dTime is zero for node1 and all other nodes that are not connected to node1
            // dTime is nonzero for all nodes that are connected to node1 (except node 1 itself)
            breadth_first_search(G,vertex(j,G), boost::visitor(vis));

            if (dtime[GetNodeID(node)] > 0)
            {
                connectedToOutputNode++;
            }

        }

    }

    return (connectedToOutputNode > 0);

}


template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SetInputNode(ChastePoint<SPATIAL_DIM> location)
{
     assert(GetNode(location)->GetNumberOfAdjoiningVessels() == 1);
     GetNode(location)->SetIsInputNode(true);
}

template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SetOutputNode(ChastePoint<SPATIAL_DIM> location)
{
    assert(GetNode(location)->GetNumberOfAdjoiningVessels() == 1);
    GetNode(location)->SetIsOutputNode(true);
}


template <unsigned SPATIAL_DIM>
void CaVascularNetwork<SPATIAL_DIM>::SaveVasculatureDataToFile(string filename)
{
    // open file to write data to
    // __________________________

    std::ofstream out(filename.c_str());

    int NumberOfPoints = 0;

    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        NumberOfPoints += GetVessel(i)->GetNumberOfSegments();
    }

    out << "# vtk DataFile Version 3.0\nvtk vasculature data\nASCII\n\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << NumberOfPoints <<" float\n";

    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
        {
            out << (GetVessel(i)->GetSegmentCoordinate(j)[0]) << " " << (GetVessel(i)->GetSegmentCoordinate(j)[1]) << " ";
            if (SPATIAL_DIM > 2)
            {
                out << (GetVessel(i)->GetSegmentCoordinate(j)[2]);
            }
            else
            {
                out << 0;
            }
            out << "\n";
        }

    }

    out << "\n\n";
    out << "LINES " << GetNumberOfVesselsInNetwork() << " " <<  GetNumberOfVesselsInNetwork() + NumberOfPoints << "\n";

    int NumberOfPointsUsed = 0;

    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetNumberOfSegments() << " ";
        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
        {
            out << NumberOfPointsUsed << " ";
            NumberOfPointsUsed++;
        }
        out << "\n";
    }


    out << "\nCELL_DATA " << GetNumberOfVesselsInNetwork() << "\n";
    out << "FIELD FieldData " << 20 + GetVessel(0)->GetNumberOfIntraVascularChemicals() << "\n";

    out << "\n";
    out << "Radius" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetRadius();
        out << "\n";
    }

    out << "\n";
    out << "UpstreamConductedStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetUpstreamConductedStimulus();
        out << "\n";
    }

    out << "\n";
    out << "DownstreamConductedStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetDownstreamConductedStimulus();
        out << "\n";
    }

    out << "\n";
    out << "ShrinkingStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetShrinkingStimulus();
        out << "\n";
    }

    out << "\n";
    out << "MetabolicStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetMetabolicStimulus();
        out << "\n";
    }

    out << "\n";
    out << "MechanicalStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetMechanicalStimulus();
        out << "\n";
    }

    out << "\n";
    out << "Viscosity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetViscosity();
        out << "\n";
    }

    out << "\n";
    out << "Impedance" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetImpedance();
        out << "\n";
    }

    out << "\n";
    out << "WallShearStress" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetWallShearStress();
        out << "\n";
    }

    out << "\n";
    out << "FlowVelocity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetFlowVelocity();
        out << "\n";
    }

    out << "\n";
    out << "FlowRate" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetFlowRate();
        out << "\n";
    }

    out << "\n";
    out << "AbsFlowVelocity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << fabs(GetVessel(i)->GetFlowVelocity());
        out << "\n";
    }

    out << "\n";
    out << "AbsFlowRate" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << fabs(GetVessel(i)->GetFlowRate());
        out << "\n";
    }

    out << "\n";
    out << "HaematocritLevel" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetHaematocritLevel();
        out << "\n";
    }

    out << "\n";
    out << "Length" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetLength();
        out << "\n";
    }

    out << "\n";
    out << "Pressure" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << 0.5*(GetVessel(i)->GetNode1()->GetPressure() + GetVessel(i)->GetNode2()->GetPressure());
        out << "\n";
    }

    out << "\n";
    out << "Tortuosity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (std::isinf(GetVessel(i)->GetTortuosity()))
        {
            // tortuosity is infinite for a circle but paraview cannot handle infinite vales so prescribe arbitrary large
            // value to print out to file if tortuosity is infinite.
            out << 10000000;
            out << "\n";
        }
        else
        {
            out << GetVessel(i)->GetTortuosity();
            out << "\n";
        }
    }

    out << "\n";
    out << "HasActiveTipCell" << " 1 " << GetNumberOfVesselsInNetwork() << " int\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->HasActiveTipCell();
        out << "\n";
    }

    out << "\n";
    out << "TimeWithLowWallShearStress" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        out << GetVessel(i)->GetTimeWithLowWallShearStress();
        out << "\n";
    }

    out << "\n";

    out << "\n";
    out << "IsPartOfNeovasculature" << " 1 " << GetNumberOfVesselsInNetwork() << " int\n";

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {

        if (GetVessel(i)->IsPartOfNeovasculature() == true)
        {
            out << "1 ";
        }
        else
        {
            out << "0 ";
        }
        out << "\n";
    }

    out << "\n";

    if (GetVessel(0)->GetNumberOfIntraVascularChemicals() > 0)
    {

        for (unsigned chemsIndex = 0; chemsIndex < GetVessel(0)->GetNumberOfIntraVascularChemicals(); chemsIndex++)
        {
            out << "\n";
            out << GetVessel(0)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetChemicalName() << "Concentration 1 " << GetNumberOfVesselsInNetwork() << " float\n";

            for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
            {

                if (GetVessel(i)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetConcentration() > 1e-15)
                {
                    out << GetVessel(i)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetConcentration();
                }
                else
                {
                    out << 0;
                }
                out << "\n";
            }



            out << "\n";
        }


    }

    out.close();

}


//template <unsigned SPATIAL_DIM>
//void CaVascularNetwork<SPATIAL_DIM>::UpdateVascularNetwork(CaBasedCellPopulation<SPATIAL_DIM>& cell_population)
//{
//
//}




// Explicit instantiation

template class CaVascularNetwork<2>;
template class CaVascularNetwork<3>;

