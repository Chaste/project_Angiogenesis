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

#include "Alarcon03HaematocritSolver.hpp"
#include "LinearSystem.hpp"
#include "VascularNode.hpp"
#include "Vessel.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::Alarcon03HaematocritSolver() : AbstractHaematocritSolver<DIM>(),
    mTHR(2.5),
    mAlpha(0.5),
    mHaematocrit(0.45)
{

}

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::~Alarcon03HaematocritSolver()
{

}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::Calculate(boost::shared_ptr<VascularNetwork<DIM> > pNetwork)
{
    // Give the vessels unique Ids
    std::vector<boost::shared_ptr<Vessel<DIM> > > vessels = pNetwork->GetVessels();
    for(unsigned idx=0; idx<vessels.size(); idx++)
    {
        vessels[idx]->SetId(idx);
    }

    // Set up the linear system
    PetscInt lhsVectorSize = vessels.size();
    unsigned max_vessels_per_branch = 5; // this could be increased if needed

    // If there are very few vessels, e.g. unit tests, set the number of non zeros to number of vessels
    if(vessels.size() < max_vessels_per_branch)
    {
        max_vessels_per_branch  = unsigned(lhsVectorSize);
    }
    LinearSystem linearSystem(lhsVectorSize, max_vessels_per_branch);

    // Currently `LinearSystem` defaults to an iterative solver for <6 dof. This check just highlights that this
    // is happening.
    if(lhsVectorSize > 6)
    {
        linearSystem.SetPcType("lu");

        // Use hypre if it is installed
        #ifndef PETSC_HAVE_HYPRE
        linearSystem.SetPcType("hypre");
        #endif //PETSC_HAVE_HYPRE
        linearSystem.SetKspType("preonly");
    }

    for(unsigned idx=0; idx<vessels.size(); idx++)
    {
        // Always have a diagonal entry for system, this sets zero haematocrit by default
        linearSystem.SetMatrixElement(idx, idx, 1);

        // Set arterial haematocrit for input nodes
        if(vessels[idx]->GetStartNode()->GetFlowProperties()->IsInputNode() or vessels[idx]->GetEndNode()->GetFlowProperties()->IsInputNode())
        {
            linearSystem.SetRhsVectorElement(idx, mHaematocrit);
        }
        // Set rhs to zero for no flow vessels. It should already be zero, but this explicitly sets it for clarity
        else if(vessels[idx]->GetFlowRate()==0.0)
        {
            linearSystem.SetRhsVectorElement(idx, 0.0);
        }
        else
        {
            // Identify the inflow node for this vessel
            boost::shared_ptr<VascularNode<DIM> > p_inflow_node;
            double flow_rate = vessels[idx]->GetFlowRate();
            if(vessels[idx]->GetFlowRate()>0)
            {
                p_inflow_node = vessels[idx]->GetStartNode();
            }
            else
            {
                p_inflow_node = vessels[idx]->GetEndNode();
            }

            // Identify the number of 'parent' and 'competitor' vessels. Parent vessels feed into the current vessel. Competitor vessels share
            // a parent vessel and do not feed in, i.e. they compete for haematocrit
            if(p_inflow_node->GetNumberOfSegments()>1)
            {
                std::vector<boost::shared_ptr<Vessel<DIM> > > parent_vessels;
                std::vector<boost::shared_ptr<Vessel<DIM> > > competitor_vessels;
                for(unsigned jdx=0; jdx<p_inflow_node->GetVesselSegments().size(); jdx++)
                {
                    // if not this vessel
                    if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()!=vessels[idx])
                    {
                        double inflow_rate = p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetFlowRate();
                        if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetEndNode()==p_inflow_node)
                        {
                            if(inflow_rate>0)
                            {
                                parent_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                            else if(inflow_rate<0)
                            {
                                competitor_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                        }
                        if(p_inflow_node->GetVesselSegment(jdx)->GetVessel()->GetStartNode()==p_inflow_node)
                        {
                            if(inflow_rate>0)
                            {
                                competitor_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                            else if(inflow_rate<0)
                            {
                                parent_vessels.push_back(p_inflow_node->GetVesselSegment(jdx)->GetVessel());
                            }
                        }
                    }
                }

                // If there are no competitor vessels the haematocrit is just the sum of the parent values
                if(competitor_vessels.size()==0)
                {
                    for(unsigned jdx=0; jdx<parent_vessels.size();jdx++)
                    {
                        linearSystem.SetMatrixElement(idx, parent_vessels[jdx]->GetId(), -1);
                    }
                }
                else
                {
                    // This is for compatibility with the paper, we could allow for more if needed
                    if(competitor_vessels.size()>1 or parent_vessels.size()>1)
                    {
                        EXCEPTION("This solver can only work with branches with connectivity 3");
                    }

                    // There is a bifurcation, apply a haematocrit splitting rule
                    double my_radius = vessels[idx]->GetRadius();
                    double competitor_radius = competitor_vessels[0]->GetRadius();
                    double my_velocity = fabs(flow_rate)/(M_PI * my_radius * my_radius);
                    double competitor_velocity = fabs(competitor_vessels[0]->GetFlowRate())/(M_PI * competitor_radius * competitor_radius);

                    if(my_velocity>mTHR*competitor_velocity)
                    {
                        // I get all the haematocrit
                        linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -1);
                    }
                    else if(my_velocity*mTHR>competitor_velocity)
                    {
                        // I get some haematocrit
                        if(my_velocity>competitor_velocity)
                        {
                            linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -1.0/(1.0+(competitor_velocity/(my_velocity*mAlpha))));
                        }
                        else if(my_velocity<competitor_velocity)
                        {
                            linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -1.0/(1.0+(mAlpha*competitor_velocity/(my_velocity))));
                        }
                        else
                        {
                            // Velocities are the same, but the partioining rule is not symmetric. This is not covered in the paper.
                            // Choose side for partitioning by vessel id so that there is some notion of conservation of H
                            if(vessels[idx]->GetId()>competitor_vessels[0]->GetId())
                            {
                                linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -1.0/(1.0+(mAlpha*competitor_velocity/(my_velocity))));
                            }
                            else
                            {
                                linearSystem.SetMatrixElement(idx, parent_vessels[0]->GetId(), -1.0/(1.0+(competitor_velocity/(my_velocity*mAlpha))));
                            }
                        }
                    }
                }
            }
        }
    }

    // Do the solve
    Vec solution = PetscTools::CreateVec(vessels.size());
    linearSystem.AssembleFinalLinearSystem();
    solution = linearSystem.Solve();
    ReplicatableVector a(solution);

    // assign haematocrit levels to vessels
    for (unsigned idx = 0; idx < vessels.size(); idx++)
    {
        for (unsigned jdx = 0; jdx < vessels[idx]->GetNumberOfSegments(); jdx++)
        {
            vessels[idx]->GetSegment(jdx)->GetFlowProperties()->SetHaematocrit(a[idx]);
        }
    }

    PetscTools::Destroy(solution);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetTHR(double THR)
{
    mTHR = THR;
    assert(mTHR > 1);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetAlpha(double Alpha)
{
    mAlpha = Alpha;
    assert(mAlpha < 1);
    assert(mAlpha > 0);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetHaematocrit(double haematocrit)
{
    mHaematocrit = haematocrit;
}

// Explicit instantiation
template class Alarcon03HaematocritSolver<2>;
template class Alarcon03HaematocritSolver<3>;
