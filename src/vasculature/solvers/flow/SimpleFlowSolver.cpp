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
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
SimpleFlowSolver<DIM>::SimpleFlowSolver()
	: mNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > >()),
	  mVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > >()),
	  mpVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> >()),
	  mNodeVesselConnectivity(std::vector<std::vector<unsigned> >()),
	  mNodeNodeConnectivity(std::vector<std::vector<unsigned> >()),
	  mBoundaryConditionNodeIndices(std::vector<unsigned> ()),
	  mUnconnectedNodeIndices(std::vector<unsigned> ()),
	  mpLinearSystem(boost::shared_ptr<LinearSystem>()),
	  mMultiplier(1.0e6),
	  mIsSetUp(false)
{

}

template<unsigned DIM>
SimpleFlowSolver<DIM>::~SimpleFlowSolver()
{

}

template<unsigned DIM>
void SimpleFlowSolver<DIM>::SetUp(boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork)
{
	mpVesselNetwork = pVascularNetwork;
	mVessels = mpVesselNetwork->GetVessels();
	mNodes = mpVesselNetwork->GetVesselEndNodes();

	unsigned num_nodes = mNodes.size();
	unsigned max_branches = mpVesselNetwork->GetMaxBranchesOnNode();

	// Set up the system
	mpLinearSystem = boost::shared_ptr<LinearSystem>(new LinearSystem (num_nodes, max_branches + 1));

	// If the network is small the preconditioner is turned off in LinearSystem,
	// so an iterative solver is used instead.
	if (num_nodes >= 6)
	{
		mpLinearSystem->SetPcType("lu");
		mpLinearSystem->SetKspType("preonly");
	}

    // Get the node-vessel and node-node connectivity
    mNodeVesselConnectivity = pVascularNetwork->GetNodeVesselConnectivity();
    mNodeNodeConnectivity = pVascularNetwork->GetNodeNodeConnectivity();

    // Get the boundary condition nodes
	std::vector<boost::shared_ptr<VascularNode<DIM> > > boundary_condition_nodes;
	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
		if (mNodes[node_index]->IsInputNode() || mNodes[node_index]->IsOutputNode())
		{
			boundary_condition_nodes.push_back(mNodes[node_index]);
			mBoundaryConditionNodeIndices.push_back(node_index);
		}
	}

	// Get the nodes that correspond to segments that are not connected to the rest of the network
    std::vector<bool> connected = mpVesselNetwork->IsConnected(boundary_condition_nodes, mNodes);
	std::vector<unsigned> mUnconnectedNodeIndices;
	for (unsigned node_index = 0; node_index < num_nodes; node_index++)
	{
	    if (!connected[node_index])
	    {
	    	mUnconnectedNodeIndices.push_back(node_index);
	    }
	}

	UpdateImpedances();

	mIsSetUp = true;
}

template<unsigned DIM>
void SimpleFlowSolver<DIM>::UpdateImpedances()
{
	mpLinearSystem->ZeroLinearSystem();

	// Get the impedances, scale them by the maximum impedance to remove small values from the system matrix
	std::vector<double> scaled_impedances;
	for (unsigned vessel_index = 0; vessel_index < mVessels.size(); vessel_index++)
	{
		double impedance = mVessels[vessel_index]->GetImpedance();
		if (impedance <= 0.0)
		{
			EXCEPTION("Impedance should be a positive number.");
		}
		scaled_impedances.push_back(impedance);
	}

    // Set up the system matrix
	for (unsigned node_index = 0; node_index < mNodes.size(); node_index++)
	{
		bool is_bc_node = (std::find(mBoundaryConditionNodeIndices.begin(), mBoundaryConditionNodeIndices.end(), node_index) != mBoundaryConditionNodeIndices.end());
		bool is_unconnected_node = (std::find(mUnconnectedNodeIndices.begin(), mUnconnectedNodeIndices.end(), node_index) != mUnconnectedNodeIndices.end());

		if(is_bc_node or is_unconnected_node)
		{
			mpLinearSystem->AddToMatrixElement(node_index, node_index, mMultiplier);
		}
		else
		{
			for (unsigned vessel_index = 0; vessel_index < mNodeVesselConnectivity[node_index].size(); vessel_index++)
			{
				double impedance = scaled_impedances[mNodeVesselConnectivity[node_index][vessel_index]];

				// Add the inverse impedances to the linear system
				mpLinearSystem->AddToMatrixElement(node_index, node_index, (-1/impedance) * mMultiplier); // Aii
				mpLinearSystem->AddToMatrixElement(node_index, mNodeNodeConnectivity[node_index][vessel_index], (+1/impedance)*mMultiplier); // Aij
			}
		}
	}

	// Update the RHS
	for (unsigned bc_index = 0; bc_index < mBoundaryConditionNodeIndices.size(); bc_index++)
	{
		mpLinearSystem->AddToRhsVectorElement(mBoundaryConditionNodeIndices[bc_index],
				mNodes[mBoundaryConditionNodeIndices[bc_index]]->GetPressure()*mMultiplier);
	}
	mpLinearSystem->AssembleIntermediateLinearSystem();
}

template<unsigned DIM>
bool SimpleFlowSolver<DIM>::IsSetUp()
{
	return mIsSetUp;
}

template<unsigned DIM>
void SimpleFlowSolver<DIM>::Implement(boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork)
{
	if (!mIsSetUp)
	{
		SetUp(pVascularNetwork);
	}

	// Assemble and solve the final system
	mpLinearSystem->AssembleFinalLinearSystem();
	Vec solution = PetscTools::CreateVec(mNodes.size());
	solution = mpLinearSystem->Solve();

	// Recover the pressure of the vessel nodes
	ReplicatableVector a(solution);
	for (unsigned node_index = 0; node_index < mNodes.size(); node_index++)
	{
		mNodes[node_index]->SetPressure(a[node_index]);
	}

	/**
		Set the segment flow rates and nodal pressures
	 */
	for (unsigned vessel_index = 0; vessel_index < mVessels.size(); vessel_index++)
	{
		double start_node_pressure = mVessels[vessel_index]->GetStartNode()->GetPressure();
		double end_node_pressure = mVessels[vessel_index]->GetEndNode()->GetPressure();
		double flow_rate = (start_node_pressure - end_node_pressure)/mVessels[vessel_index]->GetImpedance();

		std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = mVessels[vessel_index]->GetSegments();
		double pressure = start_node_pressure;
		for (unsigned segment_index = 0; segment_index < segments.size()-1; segment_index++)
		{
			pressure -= segments[segment_index]->GetImpedance() * flow_rate;
			segments[segment_index]->GetNode(1)->SetPressure(pressure);
			segments[segment_index]->SetFlowRate(flow_rate);
		}
		segments[segments.size() - 1]->SetFlowRate(flow_rate);
	}

	/*
	 * clean up
	 */
	PetscTools::Destroy(solution);
}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

