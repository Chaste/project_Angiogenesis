/*
 * SimpleFlowSolver.cpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#include "SimpleFlowSolver.hpp"
#include "PetscTools.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
SimpleFlowSolver<DIM>::SimpleFlowSolver()
{


}

template<unsigned DIM>
SimpleFlowSolver<DIM>::~SimpleFlowSolver()
{


}

template<unsigned DIM>
void SimpleFlowSolver<DIM>::Implement(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{


    /**
            @note Flow is defined as positive from node1 to node2 by convention - positive flow in a vessel
            attached to a node by its "node 0" involves flow out of the node. Conversely positive flow in
            a vessel attached to a node by its "node 1" involves flow in to the node.

            The first matrix assembled in this method assigns a coefficient to each vessel attached to a
            node which tells us whether positive flow in that vessel means flow in to or out of that node.
            If the coefficient is -1 then positive flow in the vessel means that flow is out of the node;
            if the coefficient is +1 then positive flow in the vessel means that flow is in to the node.
     */

    // assemble linear system for node pressure calculation

    PetscInt lhsVectorSize = vascularNetwork->GetNumberOfNodesInNetwork();
    unsigned rowPreallocation = 10;

    LinearSystem linearSystem(lhsVectorSize,rowPreallocation);

    // assemble structure to aid calculation

    vector< vector<int> > flowCoefficient;
    unsigned numberOfRows = vascularNetwork->GetNumberOfNodesInNetwork();
    unsigned numberOfCollumns = vascularNetwork->GetNumberOfSegmentsInNetwork();

    for (unsigned rows = 0; rows < numberOfRows; rows++)
    {
        flowCoefficient.push_back(vector<int>()); // Add an empty row
    }

    for (unsigned cols = 0; cols < numberOfCollumns; cols++)
    {
        for (unsigned rows = 0; rows < numberOfRows; rows++)
        {
            flowCoefficient[rows].push_back(0); // Add column to all rows
        }
    }

    std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetNodes();
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();
    typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator it;

    unsigned index = 0;

    for(it = nodes.begin(); it != nodes.end(); it++)
    {

        // flowCoefficient[nodeID][vesselID] = -1 if node is node 0 of vessel - denotes +ve flow if vessel flows out of node
        // flowCoefficient[nodeID][vesselID] = 1 if node is node 1 of vessel - denotes +ve flow if vessel flows into node

        for (unsigned j = 0; j < (*it)->GetNumberOfSegments(); j++)
        {
            if ((*it) == (*it)->GetVesselSegment(j)->GetNode(0) && (*it) != (*it)->GetVesselSegment(j)->GetNode(1))
            {
                flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex((*it)->GetVesselSegment(j))] = -1;
            }
            else if ((*it) == (*it)->GetVesselSegment(j)->GetNode(1) && (*it) != (*it)->GetVesselSegment(j)->GetNode(0))
            {
                flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex((*it)->GetVesselSegment(j))] = 1;
            }
            else if ((*it) == (*it)->GetVesselSegment(j)->GetNode(1) && (*it) == (*it)->GetVesselSegment(j)->GetNode(0))
            {
                // vessel loops around on itself
                flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex((*it)->GetVesselSegment(j))] = 1;
            }
            else
            {
                Exception("Node identifies vessel segment as being adjoint but node has not been correctly assigned to segment.");
            }
        }

        index++;

    }



    for (int i = 0; i < vascularNetwork->GetNumberOfNodesInNetwork(); i++)
    {
        for (int j = 0; j < vascularNetwork->GetNumberOfVesselsInNetwork(); j++)
        {
            if (CoeffMat1.x[i][j] == 1)
            {
                CoeffMat2.x[i][i] = CoeffMat2.x[i][i] - 1/(vascularNetwork->GetVessel(j)->GetImpedance());
                CoeffMat2.x[i][vascularNetwork->GetNodeID(vascularNetwork->GetVessel(j)->GetNode(0))] =  CoeffMat2.x[i][vascularNetwork->GetNodeID(vascularNetwork->GetVessel(j)->GetNode(0))] + 1/(vascularNetwork->GetVessel(j)->GetImpedance());
            }
            if (CoeffMat1.x[i][j] == -1)
            {
                CoeffMat2.x[i][i] = CoeffMat2.x[i][i] - 1/(vascularNetwork->GetVessel(j)->GetImpedance());
                CoeffMat2.x[i][vascularNetwork->GetNodeID(vascularNetwork->GetVessel(j)->GetNode(1))] = CoeffMat2.x[i][vascularNetwork->GetNodeID(vascularNetwork->GetVessel(j)->GetNode(1))] + 1/(vascularNetwork->GetVessel(j)->GetImpedance());
            }
        }
    }

    // set constant pressure coefficients in matrix for arterial input nodes and venous output nodes and assemble pressure b_vector for node pressure calculation

    for (int i = 0; i < vascularNetwork->GetNumberOfNodesInNetwork(); i++)
    {
        if (vascularNetwork->GetNode(i)->GetData<double>("Is Input"))
        {
            for (int j = 0; j < vascularNetwork->GetNumberOfNodesInNetwork(); j++)
            {
                CoeffMat2.x[i][j] = 0;
            }
            CoeffMat2.x[i][i] = 1;
            linearSystem.AddToRhsVectorElement(I,-2*pelletPermeability*VEGFConc_pellet/(Del*pelletBindingConstant));
            Pressure_bVector.x[i] = vascularNetwork->GetArterialInputPressure();
        }

        if (vascularNetwork->GetNode(i)->GetData<double>("Is Output"))
        {
            for (int j = 0; j < vascularNetwork->GetNumberOfNodesInNetwork(); j++)
            {
                CoeffMat2.x[i][j] = 0;
            }
            CoeffMat2.x[i][i] = 1;
            linearSystem.AddToRhsVectorElement(I,-2*pelletPermeability*VEGFConc_pellet/(Del*pelletBindingConstant));
            Pressure_bVector.x[i] = vascularNetwork->GetVenousOutputPressure();
        }
    }

    // set pressure

    // set pressure to zero for all nodes which are not connected to either an input node or an output node

    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;

    Graph G;

    for (int i = 0; i < vascularNetwork->GetNumberOfVesselsInNetwork(); i++)
    {
        add_edge(vascularNetwork->GetNodeID(vascularNetwork->GetVessel(i)->GetNode(0)),vascularNetwork->GetNodeID(vascularNetwork->GetVessel(i)->GetNode(1)),G);
    }

    // typedefs

    typedef boost::graph_traits<Graph>::vertices_size_type Size;

    std::vector<int> connectedToInputOrOutputNode(num_vertices(G));

    for (int j = 0; j < vascularNetwork->GetNumberOfNodesInNetwork(); j++)
    {

        if (vascularNetwork->GetNode(j)->GetData<double>("Is Input") || vascularNetwork->GetNode(j)->GetData<double>("Is Output"))
        {
            // a vector to hold the discover time property for each vertex
            std::vector<Size> dtime(num_vertices(G));

            Size time = 0;
            bfs_time_visitor<Size*>vis(&dtime[0], time);
            breadth_first_search(G,vertex(j,G), visitor(vis));

            for (int i = 0; i < vascularNetwork->GetNumberOfNodesInNetwork(); i++)
            {

                connectedToInputOrOutputNode[j] += 1;

                if (dtime[i] > 0)
                {
                    connectedToInputOrOutputNode[i] += 1;
                }

            }
        }

    }

    linearSystem.AssembleIntermediateLinearSystem();
    linearSystem.AddToMatrixElement(I,I,-2*pelletPermeability/Del);
    linearSystem.AssembleIntermediateLinearSystem();
    linearSystem.AddToMatrixElement(I,J,Do/(Del*Del));
    linearSystem.AddToRhsVectorElement(I,-2*pelletPermeability*VEGFConc_pellet/(Del*pelletBindingConstant));

    for (int i = 0; i < vascularNetwork->GetNumberOfNodesInNetwork(); i++)
    {
        if (connectedToInputOrOutputNode[i] == 0)
        {
            for (int j = 0; j < vascularNetwork->GetNumberOfNodesInNetwork(); j++)
            {
                CoeffMat2.x[i][j] = 0;
            }
            CoeffMat2.x[i][i] = 1;
        }
    }

    MatSetUp(CoeffPetsc);

    VecCreateMPI(comm,PETSC_DECIDE,matrix_size,&Rhs);
    VecDuplicate(Rhs,&Solution);

    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            if (CoeffMat2.x[i][j] != 0)
            {
                MatSetValues(CoeffPetsc,1,&i,1,&j,&CoeffMat2.x[i][j],INSERT_VALUES);
            }
        }
        VecSetValue(Rhs,i,Pressure_bVector.x[i],INSERT_VALUES);
    }

    Vec solution = PetscTools::CreateVec(vascularNetwork->GetNumberOfNodesInNetwork());
    solution = linearSystem.Solve();
    ReplicatableVector a(solution);

    // assign node pressures to vessels

    for (int i = 0; i < vascularNetwork->GetNumberOfNodesInNetwork(); i++)
    {
        vascularNetwork->GetNode(i)->SetPressure((double)a[i]);
    }

    /*
     * clean up
     */
    PetscTools::Destroy(solution);

}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

