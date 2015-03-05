/*
 * SimpleFlowSolver.cpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#include "SimpleFlowSolver.hpp"
#include "CaVesselSegment.hpp"
#include "VasculatureData.hpp"
#include "LinearSystem.hpp"
#include "Exception.hpp"
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

	PetscInt lhsVectorSize = vascularNetwork->GetNumberOfNodes();

	// assemble structure to aid calculation

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetVectorOfNodes();
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();


	std::vector< std::vector<int> > flowCoefficient;
	unsigned numberOfRows = nodes.size();
	unsigned numberOfCollumns = segments.size();

	for (unsigned rows = 0; rows < numberOfRows; rows++)
	{
		flowCoefficient.push_back(std::vector<int>()); // Add an empty row
	}

	for (unsigned cols = 0; cols < numberOfCollumns; cols++)
	{
		for (unsigned rows = 0; rows < numberOfRows; rows++)
		{
			flowCoefficient[rows].push_back(0); // Add column to all rows
		}
	}


	// calculate maximum number of non-zeros entries per row in LHS-matrix of
	// linearSystem
	unsigned max_non_zeros = 0;

	for(unsigned index = 0; index < nodes.size(); index++)
	{

		unsigned non_zeros_this_row = 1;

		// flowCoefficient[nodeID][vesselSegmentID] = -1 if node is node 0 of vessel segment - denotes +ve flow if vessel flows out of node
		// flowCoefficient[nodeID][vesselSegmentID] = 1 if node is node 1 of vessel segment - denotes +ve flow if vessel flows into node

		for (unsigned j = 0; j < nodes[index]->GetNumberOfSegments(); j++)
		{
			if (nodes[index] == nodes[index]->GetVesselSegment(j)->GetNode(0) && nodes[index] != nodes[index]->GetVesselSegment(j)->GetNode(1))
			{
				flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex(nodes[index]->GetVesselSegment(j))] = -1;
				non_zeros_this_row++;
			}
			else if (nodes[index] == nodes[index]->GetVesselSegment(j)->GetNode(1) && nodes[index] != nodes[index]->GetVesselSegment(j)->GetNode(0))
			{
				flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex(nodes[index]->GetVesselSegment(j))] = 1;
				non_zeros_this_row++;
			}
			else if (nodes[index] == nodes[index]->GetVesselSegment(j)->GetNode(1) && nodes[index] == nodes[index]->GetVesselSegment(j)->GetNode(0))
			{
				// vessel loops around on itself
				// todo check that a segment can never have two nodes that are co-located - this statement can then be deleted
				NEVER_REACHED;
				flowCoefficient[index][vascularNetwork->GetVesselSegmentIndex(nodes[index]->GetVesselSegment(j))] = 1;
			}
			else
			{
				EXCEPTION("Node identifies vessel segment as being adjoint but node has not been correctly assigned to segment.");
			}
		}

		if (non_zeros_this_row > max_non_zeros)
		{
			max_non_zeros = non_zeros_this_row;
		}

	}

	LinearSystem linearSystem(lhsVectorSize,max_non_zeros);



	for (unsigned node = 0; node < nodes.size(); node++)
	{
		for (unsigned seg = 0; seg < segments.size(); seg++)
		{
			if (flowCoefficient[node][seg] == 1)
			{
				double impedance = segments[seg]->template GetData<double>("Impedance");
				if (impedance <= 0)
				{
					EXCEPTION("Impedance should be a positive number.");
				}
				linearSystem.AddToMatrixElement(node,node,-1/impedance);
				linearSystem.AddToMatrixElement(node,vascularNetwork->GetNodeIndex(segments[seg]->GetNode(0)),+1/impedance);
			}
			if (flowCoefficient[node][seg] == -1)
			{
				double impedance = segments[seg]->template GetData<double>("Impedance");
				if (impedance <= 0)
				{
					EXCEPTION("Impedance should be a positive number.");
				}
				linearSystem.AddToMatrixElement(node,node,-1/impedance);
				linearSystem.AddToMatrixElement(node,vascularNetwork->GetNodeIndex(segments[seg]->GetNode(1)),+1/impedance);
			}
		}
	}

	linearSystem.AssembleIntermediateLinearSystem();

	// set constant pressure coefficients in matrix for arterial input nodes and venous output nodes and assemble pressure b_vector for node pressure calculation

	for (unsigned nodeIndex1 = 0; nodeIndex1 < nodes.size(); nodeIndex1++)
	{
		if (nodes[nodeIndex1]->template GetData<bool>("Is Input") || nodes[nodeIndex1]->template GetData<bool>("Is Output"))
		{
			for (unsigned nodeIndex2 = 0; nodeIndex2 < nodes.size(); nodeIndex2++)
			{
				linearSystem.SetMatrixElement(nodeIndex1,nodeIndex2,0);
			}
			linearSystem.SetMatrixElement(nodeIndex1,nodeIndex1,1);
			linearSystem.AssembleIntermediateLinearSystem();
			linearSystem.AddToRhsVectorElement(nodeIndex1,nodes[nodeIndex1]->template GetData<double>("Pressure"));
		}
	}

	// set pressure to zero for all nodes which are not connected to either an input node or an output node

	for (unsigned candidateBoundaryNode = 0; candidateBoundaryNode < nodes.size(); candidateBoundaryNode++)
	{

		if (vascularNetwork->GetNode(candidateBoundaryNode)->template GetData<bool>("Is Input") || vascularNetwork->GetNode(candidateBoundaryNode)->template GetData<bool>("Is Output"))
		{

			for (unsigned testNode = 0; testNode < vascularNetwork->GetNumberOfNodes(); testNode++)
			{

				if (!vascularNetwork->IsConnected(nodes[candidateBoundaryNode],nodes[testNode]))
				{
					for (unsigned nodeIndex1 = 0; nodeIndex1 < nodes.size(); nodeIndex1++)
					{
						linearSystem.SetMatrixElement(testNode,nodeIndex1,0);
					}
					linearSystem.SetMatrixElement(testNode,testNode,1);
				}

			}
		}

	}

	linearSystem.AssembleFinalLinearSystem();

	Vec solution = PetscTools::CreateVec(nodes.size());
	solution = linearSystem.Solve();
	ReplicatableVector a(solution);

	// assign node pressures to vessels

	for (unsigned nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++)
	{
		nodes[nodeIndex]->SetData("Pressure",a[nodeIndex]);
	}

	for (unsigned segIndex = 0; segIndex < segments.size(); segIndex++)
	{
		double flow_rate = (segments[segIndex]->GetNode(0)->template GetData<double>("Pressure")
				- segments[segIndex]->GetNode(1)->template GetData<double>("Pressure"))/segments[segIndex]->template GetData<double>("Impedance");
		segments[segIndex]->SetData("Flow Rate",flow_rate);
		segments[segIndex]->SetData("Absolute Flow Rate",fabs(flow_rate));
		segments[segIndex]->GetVessel()->SetData("Flow Rate",flow_rate);
		segments[segIndex]->GetVessel()->SetData("Absolute Flow Rate", fabs(flow_rate));
	}

	/*
	 * clean up
	 */
	PetscTools::Destroy(solution);

}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

