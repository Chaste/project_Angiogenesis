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

#include <vector>
#include <algorithm>
#include "SimpleFlowSolver.hpp"
#include "CaVesselSegment.hpp"
#include "VasculatureData.hpp"
#include "LinearSystem.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

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

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetNodes();
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
    linearSystem.SetPcType("lu");
    PetscOptionsSetValue("-pc_factor_mat_solver_package", "umfpack");
    PetscOptionsSetValue("-pc_factor_zeropivot", 0);
    linearSystem.SetKspType("gmres");

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
			double impedance = p_each_segment->GetImpedance();
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
		if (nodes[node_index]->IsInputNode() || nodes[node_index]->IsOutputNode())
		{
			rows_to_be_zerod.push_back(node_index);
			rhs_pressures.push_back(nodes[node_index]->GetPressure());
		}
	}


	linearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_be_zerod, 1.0);
	linearSystem.AssembleIntermediateLinearSystem();

	for (unsigned bc_index = 0; bc_index < rhs_pressures.size(); bc_index++)
	{
		linearSystem.AddToRhsVectorElement(rows_to_be_zerod[bc_index], rhs_pressures[bc_index]);
	}

	std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes;

    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        if (nodes[node_index]->IsInputNode() || nodes[node_index]->IsOutputNode())
        {
            sourceNodes.push_back(nodes[node_index]);
        }
     }


    std::vector<bool> connected = vascularNetwork->IsConnected(sourceNodes,nodes);

	// set pressure to zero for all nodes which are not connected to either an input node or an output node
	std::vector<unsigned> rows_to_be_zerod2;

	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
	    if (!connected[node_index])
	    {
	        rows_to_be_zerod2.push_back(node_index);
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
		nodes[node_index]->SetPressure(a[node_index]);
	}

	/**
		Set the segment flow rates
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
	{
		double start_node_pressure = segments[segment_index]->GetNode(0)->GetPressure();
		double end_node_pressure = segments[segment_index]->GetNode(1)->GetPressure();

		double flow_rate = (start_node_pressure - end_node_pressure)/segments[segment_index]->GetImpedance();
		segments[segment_index]->SetFlowRate(flow_rate);
	}

	/*
	 * clean up
	 */
	PetscTools::Destroy(solution);

}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

