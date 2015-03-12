/*
 * SimpleFlowSolver.cpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#include <vector>
#include <algorithm>
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

	/* Balance pressures in the network based on at least two defined pressures and no external flow.
	 * Through applying mass (flow) balance at nodes this results in the system: Ap = 0
	 * where pi are the nodal pressures and Aij is a matrix of inverse impedances.
	 * Aii = -sum of inverse impedances of segments attached to node i and Aij is the sum of impedances of
	 * segments between node i and node j. Practically this corresponds to only a single segment.
	 */

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetVectorOfNodes();
	unsigned num_nodes = nodes.size();

	// Get maximum number of segments attached to a node in the whole network. The system then has number of
	// segments + 1 non-zero entries.
	unsigned max_num_segments = 0;
	for(unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		boost::shared_ptr<VascularNode<DIM> > p_each_node = nodes[node_index];
		unsigned num_segments_on_node = nodes[node_index]->GetNumberOfSegments();

		if (num_segments_on_node > max_num_segments)
		{
			max_num_segments = num_segments_on_node;
		}
	}

	// Set up the system
	PetscInt lhsVectorSize = num_nodes;
	LinearSystem linearSystem(lhsVectorSize, max_num_segments + 1);

	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		boost::shared_ptr<VascularNode<DIM> > p_each_node = nodes[node_index];
		unsigned num_segments_on_node = p_each_node->GetNumberOfSegments();

		for (unsigned segment_index = 0; segment_index < num_segments_on_node; segment_index++)
		{
			boost::shared_ptr<CaVesselSegment<DIM> > p_each_segment = p_each_node->GetVesselSegment(segment_index);

			// Get the segment start and end nodes
			boost::shared_ptr<VascularNode<DIM> > p_segment_start_node = p_each_segment->GetNode(0);
			boost::shared_ptr<VascularNode<DIM> > p_segment_end_node = p_each_segment->GetNode(1);

			if(p_segment_start_node->IsCoincident(p_segment_end_node))
			{
				EXCEPTION("The network has a zero length segment. The flow problem cannot be solved. Try merging coincident nodes.");
			}

			// Get the segment impedance
			double impedance = p_each_segment->template GetData<double>("Impedance");
			if (impedance <= 0)
			{
				EXCEPTION("Impedance should be a positive number.");
			}

			unsigned node2_index;

			if (p_each_node == p_segment_start_node)
			{
				typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter = std::find(nodes.begin(), nodes.end(), p_segment_end_node);
				node2_index = std::distance(nodes.begin(), node_iter);
			}
			else
			{
				typename std::vector<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter = std::find(nodes.begin(), nodes.end(), p_segment_start_node);
				node2_index = std::distance(nodes.begin(), node_iter);
			}

			// Add the inverse impedances to the linear system
			linearSystem.AddToMatrixElement(node_index, node_index, -1/impedance); // Aii
			linearSystem.AddToMatrixElement(node_index, node2_index, +1/impedance); // Aij
		}
	}

	linearSystem.AssembleIntermediateLinearSystem();

	// set constant pressure coefficients in matrix for arterial input nodes and venous output nodes and assemble pressure b_vector for node pressure calculation
	std::vector<unsigned> rows_to_be_zerod;
	std::vector<double> rhs_pressures;
	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		if (nodes[node_index]->template GetData<bool>("Is Input") || nodes[node_index]->template GetData<bool>("Is Output"))
		{
			rows_to_be_zerod.push_back(node_index);
			rhs_pressures.push_back(nodes[node_index]->template GetData<double>("Pressure"));
		}
	}

	linearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_be_zerod, 1.0);
	linearSystem.AssembleIntermediateLinearSystem();

	for (unsigned bc_index = 0; bc_index < rhs_pressures.size(); bc_index++)
	{
		linearSystem.AddToRhsVectorElement(rows_to_be_zerod[bc_index], rhs_pressures[bc_index]);
	}

	// set pressure to zero for all nodes which are not connected to either an input node or an output node
	std::vector<unsigned> rows_to_be_zerod2;

	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		if (nodes[node_index]->template GetData<bool>("Is Input") || nodes[node_index]->template GetData<bool>("Is Output"))
		{
			for (unsigned test_node_index = 0; test_node_index < num_nodes; test_node_index++)
			{
				if (!vascularNetwork->IsConnected(nodes[node_index], nodes[test_node_index]))
				{
					rows_to_be_zerod2.push_back(test_node_index);
				}
			}
		}
	}
	linearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_be_zerod2, 1.0);

	// Assemble and solve the system
	linearSystem.AssembleFinalLinearSystem();
	Vec solution = PetscTools::CreateVec(nodes.size());
	solution = linearSystem.Solve();

	// Recover the nodal pressures
	ReplicatableVector a(solution);
	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		nodes[node_index]->SetData("Pressure", a[node_index]);
	}

	/**
		Set the segment flow rates
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
	{
		double start_node_pressure = segments[segment_index]->GetNode(0)->template GetData<double>("Pressure");
		double end_node_pressure = segments[segment_index]->GetNode(1)->template GetData<double>("Pressure");

		double flow_rate = (start_node_pressure - end_node_pressure)/segments[segment_index]->template GetData<double>("Impedance");
		segments[segment_index]->SetData("Flow Rate",flow_rate);
		segments[segment_index]->SetData("Absolute Flow Rate",fabs(flow_rate));
		segments[segment_index]->GetVessel()->SetData("Flow Rate",flow_rate);
		segments[segment_index]->GetVessel()->SetData("Absolute Flow Rate", fabs(flow_rate));
	}

	/*
	 * clean up
	 */
	PetscTools::Destroy(solution);

}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

