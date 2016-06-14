///*
//
//Copyright (c) 2005-2015, University of Oxford.
//All rights reserved.
//
//University of Oxford means the Chancellor, Masters and Scholars of the
//University of Oxford, having an administrative office at Wellington
//Square, Oxford OX1 2JD, UK.
//
//This file is part of Chaste.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
//   contributors may be used to endorse or promote products derived from this
//   software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS unsignedERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "BetteridgeHaematocritSolver.hpp"
#include "LinearSystem.hpp"
#include "VascularNode.hpp"
#include "Vessel.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
BetteridgeHaematocritSolver<DIM>::BetteridgeHaematocritSolver() : AbstractHaematocritSolver<DIM>(),
    mTHR(2.5),
    mAlpha(0.5),
    mHaematocrit(0.45)
{

}

template<unsigned DIM>
BetteridgeHaematocritSolver<DIM>::~BetteridgeHaematocritSolver()
{

}

template<unsigned DIM>
void BetteridgeHaematocritSolver<DIM>::SetTHR(units::quantity<unit::dimensionless> THR)
{
    mTHR = THR;
}

template<unsigned DIM>
void BetteridgeHaematocritSolver<DIM>::SetAlpha(units::quantity<unit::dimensionless> Alpha)
{
    mAlpha = Alpha;
}

template<unsigned DIM>
void BetteridgeHaematocritSolver<DIM>::SetHaematocrit(units::quantity<unit::dimensionless> haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void BetteridgeHaematocritSolver<DIM>::Calculate(boost::shared_ptr<VascularNetwork<DIM> > pNetwork)
{
    // Give the vessels unique Ids
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels = pNetwork->GetVessels();
    for(unsigned idx=0; idx<vessels.size(); idx++)
    {
        vessels[idx]->SetId(idx);
    }

    // Set up the linear system
    PetscInt lhsVectorSize = vessels.size();
    unsigned max_vessels_per_branch = 5;
    if(vessels.size() < max_vessels_per_branch)
    {
        max_vessels_per_branch  = unsigned(lhsVectorSize);
    }
    LinearSystem linearSystem(lhsVectorSize, max_vessels_per_branch);
    if(lhsVectorSize > 6)
    {
        linearSystem.SetPcType("lu");
        #ifndef PETSC_HAVE_HYPRE
        linearSystem.SetPcType("hypre");
        #endif //PETSC_HAVE_HYPRE
        linearSystem.SetKspType("preonly");
    }

    std::vector<std::vector<unsigned> > update_indices;
    for(unsigned idx=0; idx<vessels.size(); idx++)
    {
        // Always have a diagonal entry for system, this sets zero haematocrit by default
        linearSystem.SetMatrixElement(idx, idx, 1);
        if(vessels[idx]->GetStartNode()->GetFlowProperties()->IsInputNode() or vessels[idx]->GetEndNode()->GetFlowProperties()->IsInputNode())
        {
            linearSystem.SetRhsVectorElement(idx, mHaematocrit);
        }
        // Set rhs to zero, it should already be zero but this explicitly captures the no flow case
        else if(vessels[idx]->GetFlowRate()==0.0*unit::unit_flow_rate)
        {
            linearSystem.SetRhsVectorElement(idx, 0.0);
        }
        else
        {
            // Identify inflow node
            boost::shared_ptr<VascularNode<DIM> > p_inflow_node;
            units::quantity<unit::flow_rate> flow_rate= vessels[idx]->GetFlowRate();
            if(vessels[idx]->GetFlowRate()>0 * unit::unit_flow_rate)
            {
                p_inflow_node = vessels[idx]->GetStartNode();
            }
            else
            {
                p_inflow_node = vessels[idx]->GetEndNode();
            }

            // Identify number of inflow and outflow vessels
            if(p_inflow_node->GetNumberOfSegments()>1)
            {
                std::vector<boost::shared_ptr<Vessel<DIM> > > parent_vessels;
                std::vector<boost::shared_ptr<Vessel<DIM> > > competitor_vessels;
                for(unsigned jdx=0; jdx<p_inflow_node->GetVesselSegments().size(); jdx++)
                {
                    // if not this vessel
                    if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()!=vessels[idx])
                    {
                        units::quantity<unit::flow_rate> inflow_rate = p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetFlowRate();
                        if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetEndNode()==p_inflow_node)
                        {
                            if(inflow_rate>0.0 * unit::unit_flow_rate)
                            {
                                parent_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                            else if(inflow_rate<0.0 * unit::unit_flow_rate)
                            {
                                competitor_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                        }
                        if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetStartNode()==p_inflow_node)
                        {
                            if(inflow_rate>0.0 * unit::unit_flow_rate)
                            {
                                competitor_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                            else if(inflow_rate<0.0 * unit::unit_flow_rate)
                            {
                                parent_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                        }
                    }
                }

                // If there are no competitor vessels the haematocrit is just the sum of the parent values
                if(competitor_vessels.size()==0 or units::fabs(competitor_vessels[0]->GetFlowRate()) == 0.0 * unit::unit_flow_rate)
                {
                    for(unsigned jdx=0; jdx<parent_vessels.size();jdx++)
                    {
                        linearSystem.SetMatrixElement(idx, parent_vessels[jdx]->GetId(), -fabs(parent_vessels[jdx]->GetFlowRate()/flow_rate));
                    }
                }
                else
                {
                    if(competitor_vessels.size()>1 or parent_vessels.size()>1)
                    {
                        EXCEPTION("This solver can only work with branches with connectivity 3");
                    }

                    // There is a bifurcation, apply a haematocrit splitting rule
                    units::quantity<unit::length> my_radius = vessels[idx]->GetRadius();
                    units::quantity<unit::length> competitor_radius = competitor_vessels[0]->GetRadius();
                    units::quantity<unit::velocity> my_velocity = units::fabs(flow_rate)/(M_PI * my_radius * my_radius);
                    units::quantity<unit::velocity> competitor_velocity = fabs(competitor_vessels[0]->GetFlowRate())/(M_PI * competitor_radius * competitor_radius);

                    // Alpha now depends on haematocrit in the parent vessel, so system is non-linear. Will solve iteratively
                    units::quantity<unit::dimensionless> alpha = 1.0 - parent_vessels[0]->GetHaematocrit();
                    units::quantity<unit::dimensionless> term = alpha * (my_velocity/competitor_velocity-1.0);
                    units::quantity<unit::dimensionless> flow_ratio_pm = fabs(parent_vessels[0]->GetFlowRate())/fabs(flow_rate);
                    units::quantity<unit::dimensionless> flow_ratio_cm = fabs(competitor_vessels[0]->GetFlowRate())/fabs(flow_rate);

                    double numer = flow_ratio_pm;
                    double denom = 1.0+flow_ratio_cm*(1.0/(1.0+term));
                    linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -numer/denom);

                    // Save the indices for later updating
                    std::vector<unsigned> local_update_indics = std::vector<unsigned>(3);
                    local_update_indics[0] = idx;
                    local_update_indics[1] = parent_vessels[0]->GetId();
                    local_update_indics[2] = competitor_vessels[0]->GetId();
                    update_indices.push_back(local_update_indics);
                }
            }
        }
    }

    // Set the parameters for iteration
    double tolerance = 1e-3;
    double max_iterations = 1000;
    double residual = DBL_MAX;
    int iterations = 0;

    while(residual > tolerance && iterations < max_iterations)
    {
        if(iterations>0 and update_indices.size()>0)
        {
            // Update the system
            linearSystem.SwitchWriteModeLhsMatrix();
            for(unsigned idx=0; idx<update_indices.size();idx++)
            {
                units::quantity<unit::flow_rate> flow_rate = vessels[update_indices[idx][0]]->GetFlowRate();
                units::quantity<unit::length> my_radius = vessels[update_indices[idx][0]]->GetRadius();
                units::quantity<unit::length> competitor_radius = vessels[update_indices[idx][2]]->GetRadius();
                units::quantity<unit::velocity> my_velocity = fabs(flow_rate)/(M_PI * my_radius * my_radius);
                units::quantity<unit::velocity> competitor_velocity = fabs(vessels[update_indices[idx][2]]->GetFlowRate())/(M_PI * competitor_radius * competitor_radius);
                units::quantity<unit::dimensionless> alpha = 1.0 - vessels[update_indices[idx][1]]->GetHaematocrit();

                double term = alpha * (my_velocity/competitor_velocity-1.0);
                double flow_ratio_pm = fabs(vessels[update_indices[idx][1]]->GetFlowRate()/flow_rate);
                double flow_ratio_cm = fabs(vessels[update_indices[idx][2]]->GetFlowRate()/flow_rate);
                double numer = flow_ratio_pm;
                double denom = 1.0+flow_ratio_cm*(1.0/(1.0+term));

                linearSystem.SetMatrixElement(update_indices[idx][0], update_indices[idx][1], -numer/denom);
            }
        }

        Vec solution = PetscTools::CreateVec(vessels.size());
        linearSystem.AssembleFinalLinearSystem();
        solution = linearSystem.Solve();
        ReplicatableVector a(solution);

        // Get the residual
        residual = 0;
        for (unsigned i = 0; i < vessels.size(); i++)
        {
            if(vessels[i]->GetHaematocrit() - a[i] > residual)
            {
                residual = vessels[i]->GetHaematocrit() - a[i];
            }
        }

        // assign haematocrit levels to vessels
        for (unsigned idx = 0; idx < vessels.size(); idx++)
        {
            for (unsigned jdx = 0; jdx < vessels[idx]->GetNumberOfSegments(); jdx++)
            {
                vessels[idx]->GetSegment(jdx)->GetFlowProperties()->SetHaematocrit(a[idx]);
            }
        }

        iterations++;
        if(iterations == max_iterations)
        {
            EXCEPTION("Haematocrit calculation failed to converge.");
        }

        PetscTools::Destroy(solution);
    }
}
// Explicit instantiation
template class BetteridgeHaematocritSolver<2>;
template class BetteridgeHaematocritSolver<3>;
