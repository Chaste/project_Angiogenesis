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
	: mVesselArray(),
	  mNodeArray(),
	  mArterialHaematocritLevel(0.45),
	  mArterialInputPressure(27*(1.01*pow(10.0,5)/760)),
	  mVenousOutputPressure(15*(1.01*pow(10.0,5)/760))
{
}

template <unsigned DIM>
CaVascularNetwork<DIM>::~CaVascularNetwork()
{
}

template <unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > CaVascularNetwork<DIM>::shared()
{
    return this->shared_from_this();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetVesselID(boost::shared_ptr<CaVessel<DIM> > vessel)
{
    bool vessel_found = false;
    unsigned i = 0;

    for (i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        if (vessel == mVesselArray[i])
        {
        	vessel_found = true;
            break;
        }
    }

    assert(vessel_found);
    return i;
}


template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNodeID(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{
    bool node_found = false;
    unsigned i = 0;

    for (i = 0; i < GetNumberOfNodesInNetwork(); i++)
    {
        if (node == mNodeArray[i])
        {
        	node_found = true;
            break;
        }
    }

    assert(node_found);
    return i;
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfVesselsInNetwork()
{
    return mVesselArray.size();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfNodesInNetwork()
{
    return mNodeArray.size();
}

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::GetNumberOfVesselsAtLocation(ChastePoint<DIM> coord)
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

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetVessel(int vessel_id)
{
    return mVesselArray[vessel_id];
}

template <unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVascularNetwork<DIM>::GetVessel(ChastePoint<DIM> coord, int positionInContainer)
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

template <unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVascularNetwork<DIM>::GetNode(int node_id)
{
    return mNodeArray[node_id];
}

template <unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVascularNetwork<DIM>::GetNode(ChastePoint<DIM> location)
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

template <unsigned DIM>
double CaVascularNetwork<DIM>::GetArterialHaematocritLevel()
{
    return mArterialHaematocritLevel;
}

template <unsigned DIM>
double CaVascularNetwork<DIM>::GetArterialInputPressure()
{
    return mArterialInputPressure;
}

template <unsigned DIM>
double CaVascularNetwork<DIM>::GetVenousOutputPressure()
{
    return mVenousOutputPressure;
}

//template <unsigned DIM>
//double CaVascularNetwork<DIM>::GetMeanVesselLengthOfNeovasculature()
//{
//    double totalVesselLength = 0;
//    int numberofNewVessels = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (GetVessel(i)->IsPartOfNeovasculature())
//        {
//            totalVesselLength += GetVessel(i)->GetLength();
//            numberofNewVessels++;
//        }
//
//    }
//
//    return (totalVesselLength/(double)numberofNewVessels);
//}
//
//template <unsigned DIM>
//int CaVascularNetwork<DIM>::GetNumberOfVesselsByLength(double lowerBoundLength, double upperBoundLength)
//{
//    int numberOfVesselsInsideRange = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (GetVessel(i)->GetLength() > lowerBoundLength && GetVessel(i)->GetLength() <= upperBoundLength)
//        {
//            numberOfVesselsInsideRange++;
//        }
//    }
//
//    return numberOfVesselsInsideRange;
//}
//
//template <unsigned DIM>
//int CaVascularNetwork<DIM>::GetNumberOfVesselsByRadius(double lowerBoundRadius, double upperBoundRadius)
//{
//    int numberOfVesselsInsideRange = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (GetVessel(i)->GetRadius() > lowerBoundRadius && GetVessel(i)->GetRadius() <= upperBoundRadius)
//        {
//            numberOfVesselsInsideRange++;
//        }
//    }
//
//    return numberOfVesselsInsideRange;
//}
//
//template <unsigned DIM>
//int CaVascularNetwork<DIM>::GetNumberOfVesselsByTortuosity(double lowerBoundTortuosity, double upperBoundTortuosity)
//{
//    int numberOfVesselsInsideRange = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (GetVessel(i)->GetTortuosity() > lowerBoundTortuosity && GetVessel(i)->GetTortuosity() <= upperBoundTortuosity)
//        {
//            numberOfVesselsInsideRange++;
//        }
//    }
//
//    return numberOfVesselsInsideRange;
//}
//
//template <unsigned DIM>
//double CaVascularNetwork<DIM>::GetMeanVesselRadiusOfNeovasculature()
//{
//    double totalVesselRadius = 0;
//    int numberofNewVessels = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (GetVessel(i)->IsPartOfNeovasculature())
//        {
//            totalVesselRadius += GetVessel(i)->GetRadius();
//            numberofNewVessels++;
//        }
//
//    }
//
//    return (totalVesselRadius/(double)numberofNewVessels);
//}
//
//template <unsigned DIM>
//double CaVascularNetwork<DIM>::GetMeanVesselTortuosityOfNeovasculature()
//{
//    double totalVesselTortuosity = 0;
//    int numberofNewVessels = 0;
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (std::isinf(GetVessel(i)->GetTortuosity()))
//        {
//            // we ignore vessels with infinite tortuosity, i.e. self-loops
//        }
//        else
//        {
//            if (GetVessel(i)->IsPartOfNeovasculature())
//            {
//                totalVesselTortuosity += GetVessel(i)->GetTortuosity();
//                numberofNewVessels++;
//            }
//        }
//    }
//
//    return (totalVesselTortuosity/(double)numberofNewVessels);
//}

template <unsigned DIM>
std::vector<boost::shared_ptr<CaVessel<DIM> > > CaVascularNetwork<DIM>::GetVessels()
{
    return mVesselArray;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetArterialHaematocritLevel(double value)
{
    mArterialHaematocritLevel = value;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetArterialInputPressure(double value)
{
    mArterialInputPressure = value;
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetVenousOutputPressure(double value)
{
    mVenousOutputPressure = value;
}

template <unsigned DIM>
bool CaVascularNetwork<DIM>::NodePresentAtLocation(ChastePoint<DIM> location)
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

template <unsigned DIM>
unsigned CaVascularNetwork<DIM>::NumberOfNodesPresentAtLocation(ChastePoint<DIM> location)
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

template <unsigned DIM>
void CaVascularNetwork<DIM>::AddVessel(boost::shared_ptr<CaVessel<DIM> > vessel)
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

template <unsigned DIM>
bool CaVascularNetwork<DIM>::VesselIsInNetwork(boost::shared_ptr<CaVessel<DIM> > vessel)
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

template <unsigned DIM>
bool CaVascularNetwork<DIM>::NodeIsInNetwork(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
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


///\ todo this could be made more general by passing in vectors of source nodes and target nodes and returning true for any targets connected
// to a source. Would avoid graph reconstruction for each query and results in only one overall method.
template <unsigned DIM>
bool CaVascularNetwork<DIM>::Connected(boost::shared_ptr<CaVascularNetworkNode<DIM> > node1, boost::shared_ptr<CaVascularNetworkNode<DIM> > node2)
{

    if (node1 == node2)
    {
        return true;
    }

    // construct graph representation of vessel network
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(GetNodeID(GetVessel(i)->GetNode1()), GetNodeID(GetVessel(i)->GetNode2()), G);
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

template <unsigned DIM>
bool CaVascularNetwork<DIM>::ConnectedToInputNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{

    if (node->GetBooleanData("IsInputNode"))
    {
        return true;
    }

    // construct graph representation of vessel network
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;

    Graph G;

    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(GetNodeID(GetVessel(i)->GetNode1()),GetNodeID(GetVessel(i)->GetNode2()), G);
    }

    // typedefs
    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    int connectedToInputNode = 0;

    for (unsigned j = 0; j < GetNumberOfNodesInNetwork(); j++)
    {

        if (GetNode(j)->GetBooleanData("IsInputNode"))
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


template <unsigned DIM>
bool CaVascularNetwork<DIM>::ConnectedToOutputNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{

    if (node->GetBooleanData("IsOutputNode"))
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

        if (GetNode(j)->GetBooleanData("IsOutputNode"))
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

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetInputNode(ChastePoint<DIM> location)
{
     assert(GetNode(location)->GetNumberOfAdjoiningVessels() == 1);
     GetNode(location)->SetBooleanData("IsInputNode", true);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SetOutputNode(ChastePoint<DIM> location)
{
    assert(GetNode(location)->GetNumberOfAdjoiningVessels() == 1);
    GetNode(location)->SetBooleanData("IsOutputNode", true);
}

template <unsigned DIM>
void CaVascularNetwork<DIM>::SaveVasculatureDataToFile(string filename)
{
    // open file to write data to
    // __________________________

	///\ todo replace with vtk polydata writer

//    std::ofstream out(filename.c_str());
//
//    int NumberOfPoints = 0;
//
//    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        NumberOfPoints += GetVessel(i)->GetNumberOfSegments();
//    }
//
//    out << "# vtk DataFile Version 3.0\nvtk vasculature data\nASCII\n\n";
//    out << "DATASET POLYDATA\n";
//    out << "POINTS " << NumberOfPoints <<" float\n";
//
//    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
//        {
//            out << (GetVessel(i)->GetSegmentCoordinate(j)[0]) << " " << (GetVessel(i)->GetSegmentCoordinate(j)[1]) << " ";
//            if (DIM > 2)
//            {
//                out << (GetVessel(i)->GetSegmentCoordinate(j)[2]);
//            }
//            else
//            {
//                out << 0;
//            }
//            out << "\n";
//        }
//
//    }
//
//    out << "\n\n";
//    out << "LINES " << GetNumberOfVesselsInNetwork() << " " <<  GetNumberOfVesselsInNetwork() + NumberOfPoints << "\n";
//
//    int NumberOfPointsUsed = 0;
//
//    for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetNumberOfSegments() << " ";
//        for (unsigned j = 0; j < GetVessel(i)->GetNumberOfSegments(); j++)
//        {
//            out << NumberOfPointsUsed << " ";
//            NumberOfPointsUsed++;
//        }
//        out << "\n";
//    }
//
//
//    out << "\nCELL_DATA " << GetNumberOfVesselsInNetwork() << "\n";
//    out << "FIELD FieldData " << 20 + GetVessel(0)->GetNumberOfIntraVascularChemicals() << "\n";
//
//    out << "\n";
//    out << "Radius" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetRadius();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "UpstreamConductedStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetUpstreamConductedStimulus();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "DownstreamConductedStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetDownstreamConductedStimulus();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "ShrinkingStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetShrinkingStimulus();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "MetabolicStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetMetabolicStimulus();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "MechanicalStimulus" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetMechanicalStimulus();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "Viscosity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetViscosity();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "Impedance" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetImpedance();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "WallShearStress" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetWallShearStress();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "FlowVelocity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetFlowVelocity();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "FlowRate" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetFlowRate();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "AbsFlowVelocity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << fabs(GetVessel(i)->GetFlowVelocity());
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "AbsFlowRate" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << fabs(GetVessel(i)->GetFlowRate());
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "HaematocritLevel" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetHaematocritLevel();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "Length" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetLength();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "Pressure" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << 0.5*(GetVessel(i)->GetNode1()->GetPressure() + GetVessel(i)->GetNode2()->GetPressure());
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "Tortuosity" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        if (std::isinf(GetVessel(i)->GetTortuosity()))
//        {
//            // tortuosity is infinite for a circle but paraview cannot handle infinite vales so prescribe arbitrary large
//            // value to print out to file if tortuosity is infinite.
//            out << 10000000;
//            out << "\n";
//        }
//        else
//        {
//            out << GetVessel(i)->GetTortuosity();
//            out << "\n";
//        }
//    }
//
//    out << "\n";
//    out << "HasActiveTipCell" << " 1 " << GetNumberOfVesselsInNetwork() << " int\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->HasActiveTipCell();
//        out << "\n";
//    }
//
//    out << "\n";
//    out << "TimeWithLowWallShearStress" << " 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//        out << GetVessel(i)->GetTimeWithLowWallShearStress();
//        out << "\n";
//    }
//
//    out << "\n";
//
//    out << "\n";
//    out << "IsPartOfNeovasculature" << " 1 " << GetNumberOfVesselsInNetwork() << " int\n";
//
//    for (unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//    {
//
//        if (GetVessel(i)->IsPartOfNeovasculature() == true)
//        {
//            out << "1 ";
//        }
//        else
//        {
//            out << "0 ";
//        }
//        out << "\n";
//    }
//
//    out << "\n";
//
//    if (GetVessel(0)->GetNumberOfIntraVascularChemicals() > 0)
//    {
//
//        for (unsigned chemsIndex = 0; chemsIndex < GetVessel(0)->GetNumberOfIntraVascularChemicals(); chemsIndex++)
//        {
//            out << "\n";
//            out << GetVessel(0)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetChemicalName() << "Concentration 1 " << GetNumberOfVesselsInNetwork() << " float\n";
//
//            for(unsigned i = 0; i < GetNumberOfVesselsInNetwork(); i++)
//            {
//
//                if (GetVessel(i)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetConcentration() > 1e-15)
//                {
//                    out << GetVessel(i)->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[chemsIndex].GetConcentration();
//                }
//                else
//                {
//                    out << 0;
//                }
//                out << "\n";
//            }
//
//
//
//            out << "\n";
//        }
//
//
//    }
//
//    out.close();
}

// Explicit instantiation
template class CaVascularNetwork<2>;
template class CaVascularNetwork<3>;

