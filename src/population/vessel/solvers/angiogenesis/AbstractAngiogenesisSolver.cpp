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

#include <boost/lexical_cast.hpp>
#include "UblasVectorInclude.hpp"
#include "UblasIncludes.hpp"
#include "RandomNumberGenerator.hpp"
#include "VascularNode.hpp"
#include "SimpleFlowSolver.hpp"
#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#include "AbstractAngiogenesisSolver.hpp"

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::AbstractAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork) :
        mpNetwork(pNetwork),
        mGrowthVelocity(10.0),
        mEndTime(10.0),
        mOutputFrequency(1),
        mOutputDirectory(),
        mNodeAnastamosisRadius(0.0),
        mPdeSolvers(),
        mGrowthDirectionModifiers(),
        mpSproutingRule(),
        mpFlowSolver(),
        mpStructuralAdaptationSolver()
{

}

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::~AbstractAngiogenesisSolver()
{

}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::AddGrowthDirectionModifier(boost::shared_ptr<AbstractGrowthDirectionModifier<DIM> > pModifier)
{
    mGrowthDirectionModifiers.push_back(pModifier);
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::AddPdeSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pPdeSolver)
{
    mPdeSolvers.push_back(pPdeSolver);
}

template<unsigned DIM>
std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > AbstractAngiogenesisSolver<DIM>::GetPdeSolvers()
{
    return mPdeSolvers;
}

template<unsigned DIM>
c_vector<double, DIM> AbstractAngiogenesisSolver<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                                                          boost::shared_ptr<VascularNode<DIM> > pNode)
{
    c_vector<double,DIM> new_direction = currentDirection;

    // Loop through the growth direction modifiers and add up the contributions
    for(unsigned idx=0; idx<mGrowthDirectionModifiers.size(); idx++)
    {
        new_direction += mGrowthDirectionModifiers[idx]->GetStrength() * mGrowthDirectionModifiers[idx]->GetGrowthDirection(currentDirection, pNode);
        new_direction /= norm_2(new_direction);
    }

    return new_direction;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetEndTime(double time)
{
    mEndTime = time;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetFlowSolver(boost::shared_ptr<SimpleFlowSolver<DIM> > pFlowSolver)
{
    mpFlowSolver = pFlowSolver;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetGrowthVelocity(double velocity)
{
    mGrowthVelocity = velocity;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetOutputDirectory(const std::string& rDirectory)
{
    mOutputDirectory = rDirectory;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetOutputFrequency(unsigned frequency)
{
    mOutputFrequency = frequency;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetSproutingRule(boost::shared_ptr<AbstractSproutingRule<DIM> > pSproutingRule)
{
    mpSproutingRule = pSproutingRule;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetStructuralAdaptationSolver(boost::shared_ptr<SimpleStructuralAdaptationSolver<DIM> > pStructuralAdaptationSolver)
{
    mpStructuralAdaptationSolver = pStructuralAdaptationSolver;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::DoSprouting()
{
    mpNetwork->UpdateNodes();

    // Get the candidate nodes
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    // Get the sprout indices and directions
    mpSproutingRule->SetNodes(nodes);
    std::vector<bool> sprout_indices = mpSproutingRule->WillSprout();
    std::vector<c_vector<double, DIM> > sprout_directions = mpSproutingRule->GetSproutDirection();

    // Do the sprouting
    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        if(sprout_indices[idx])
        {
            // Ensure it is not along the segment vectors
            bool is_along_segment_1 = std::abs(inner_prod(sprout_directions[idx],nodes[idx]->GetVesselSegments()[0]->GetUnitTangent())/
                    (norm_2(sprout_directions[idx])*norm_2(nodes[idx]->GetVesselSegments()[0]->GetUnitTangent()))) > 1 - 1.e-6;
            bool is_along_segment_2 = std::abs(inner_prod(sprout_directions[idx],nodes[idx]->GetVesselSegments()[1]->GetUnitTangent())/
                    (norm_2(sprout_directions[idx])*norm_2(nodes[idx]->GetVesselSegments()[1]->GetUnitTangent()))) > 1 - 1.e-6;

            if(!is_along_segment_1 && !is_along_segment_2)
            {
                mpNetwork->FormSprout(nodes[idx]->GetLocation(), ChastePoint<DIM>(nodes[idx]->GetLocationVector() +
                                                                                  mGrowthVelocity*sprout_directions[idx]));
            }
        }
    }
    mpNetwork->UpdateSegments();
    mpNetwork->UpdateNodes();
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::UpdateNodalPositions(const std::string& speciesLabel)
{
    mpNetwork->UpdateNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments()==1)
        {
            // Get the segment direction vector
            c_vector<double,DIM> direction = nodes[idx]->GetLocationVector() -
                 nodes[idx]->GetVesselSegment(0)->GetOppositeNode(nodes[idx])->GetLocationVector();
            direction /= norm_2(direction);

            // Create a new segment along the growth vector
            boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(nodes[idx]);
            p_new_node->SetLocation(nodes[idx]->GetLocationVector() + mGrowthVelocity * GetGrowthDirection(direction, nodes[idx]));

            if(nodes[idx]->GetVesselSegment(0)->GetVessel()->GetStartNode() == nodes[idx])
            {
                boost::shared_ptr<CaVesselSegment<DIM> > p_segment = CaVesselSegment<DIM>::Create(p_new_node, nodes[idx]);
                p_segment->SetFlowProperties(*nodes[idx]->GetVesselSegment(0)->GetFlowProperties());
                p_segment->SetRadius(nodes[idx]->GetVesselSegment(0)->GetRadius());
                nodes[idx]->GetVesselSegment(0)->GetVessel()->AddSegment(p_segment);
            }
            else
            {
                boost::shared_ptr<CaVesselSegment<DIM> > p_segment = CaVesselSegment<DIM>::Create(nodes[idx], p_new_node);
                p_segment->SetFlowProperties(*nodes[idx]->GetVesselSegment(0)->GetFlowProperties());
                p_segment->SetRadius(nodes[idx]->GetVesselSegment(0)->GetRadius());
                nodes[idx]->GetVesselSegment(0)->GetVessel()->AddSegment(p_segment);
            }
            nodes[idx]->SetIsMigrating(false);
            p_new_node->SetIsMigrating(true);
        }
    }
    mpNetwork->UpdateSegments();
    mpNetwork->UpdateNodes();
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::DoAnastamosis()
{
    // Do tip-tip anastamosis and tip-stalk anastamosis for nearby nodes
    mpNetwork->UpdateNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > moved_nodes = mpNetwork->GetNodes();
    for(unsigned idx = 0; idx < moved_nodes.size(); idx++)
    {
        if(moved_nodes[idx]->IsMigrating() && moved_nodes[idx]->GetNumberOfSegments()==1)
        {
            std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_pair = mpNetwork->GetNearestSegment(moved_nodes[idx]);
            if(segment_pair.second <= mNodeAnastamosisRadius)
            {
                // Divide the parent vessel if neccessary and set all involved nodes to non-migrating
                boost::shared_ptr<VascularNode<DIM> > p_merge_node = mpNetwork->DivideVessel(segment_pair.first->GetVessel(), moved_nodes[idx]->GetLocation());
                p_merge_node->SetIsMigrating(false);
                moved_nodes[idx]->SetIsMigrating(false);
            }
        }
    }

    // Check for crossing segments (should also do overlapping ones)
    mpNetwork->MergeCoincidentNodes();
    mpNetwork->UpdateNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > remaining_nodes = mpNetwork->GetNodes();
    for(unsigned idx = 0; idx < remaining_nodes.size(); idx++)
    {
        if(remaining_nodes[idx]->IsMigrating() && remaining_nodes[idx]->GetNumberOfSegments()==1)
        {
            std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_pair = mpNetwork->GetNearestSegment(remaining_nodes[idx]->GetVesselSegment(0));
            if(segment_pair.second <= mNodeAnastamosisRadius)
            {
                c_vector<double, DIM> divide_location = segment_pair.first->GetPointProjection(remaining_nodes[idx]->GetLocation());
                boost::shared_ptr<VascularNode<DIM> > p_merge_node =
                        mpNetwork->DivideVessel(segment_pair.first->GetVessel(), divide_location);
                p_merge_node->SetIsMigrating(false);
                // todo need to remove the overlapping segment here, otherwise will have zero length.
                remaining_nodes[idx]->SetLocation(divide_location);
                remaining_nodes[idx]->SetIsMigrating(false);
            }
        }
    }
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::Increment()
{
    unsigned num_steps = SimulationTime::Instance()->GetTimeStepsElapsed();

    mpNetwork->MergeCoincidentNodes();
    mpNetwork->UpdateVesselIds();
    mpNetwork->UpdateVesselNodes();

    // If there is a structural adaptation or flow problem solve it
    if(mpStructuralAdaptationSolver)
    {
            EXCEPTION("Structural Adapation is not implemented in this solver yet.");
    }
    else if(mpFlowSolver)
    {
        PoiseuilleImpedanceCalculator<DIM> impedance_calculator;
        impedance_calculator.Calculate(mpNetwork);
        mpFlowSolver->SetUp(mpNetwork);
        mpFlowSolver->Implement(mpNetwork);
    }

    // If there are PDEs solve them
    if(mPdeSolvers.size()>0)
    {
        for(unsigned idx=0; idx<mPdeSolvers.size(); idx++)
        {
            mPdeSolvers[idx]->SetVesselNetwork(mpNetwork);
            mPdeSolvers[idx]->SetWorkingDirectory(mOutputDirectory);
            std::string species_name = mPdeSolvers[idx]->GetPde()->GetVariableName();

            mPdeSolvers[idx]->SetFileName("/" + species_name +"_solution_" + boost::lexical_cast<std::string>(num_steps)+".vti");

            // Take the previous pde solution if needed
            if(idx>0)
            {
                for(unsigned jdx=0; jdx<mPdeSolvers[idx]->GetPde()->GetDiscreteSources().size(); jdx++)
                {
                    if(mPdeSolvers[idx]->GetPde()->GetDiscreteSources()[jdx]->GetType()==SourceType::SOLUTION)
                    {
                        mPdeSolvers[idx]->GetPde()->GetDiscreteSources()[jdx]->SetSolution(mPdeSolvers[idx-1]->GetSolution());
                    }
                }
            }

            if(mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
            {
                mPdeSolvers[idx]->Solve(true);
            }
            else
            {
                mPdeSolvers[idx]->Solve(false);
            }
        }
    }

    // Move any migrating nodes
    UpdateNodalPositions();
    DoAnastamosis();

    if(mpSproutingRule)
    {
        DoSprouting();
    }

    // Do anastamosis
    DoAnastamosis();
    mpNetwork->MergeCoincidentNodes();

    if(mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
    {
        mpNetwork->UpdateVesselIds();
        mpNetwork->Write(mOutputDirectory + "/VesselNetwork_inc_" + boost::lexical_cast<std::string>(num_steps)+".vtp");
    }
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::Run()
{
    mpNetwork->MergeCoincidentNodes();
    mpNetwork->UpdateVesselIds();
    mpNetwork->UpdateVesselNodes();
    mpNetwork->Write(mOutputDirectory + "/InputVesselNetwork.vtp");

    while(SimulationTime::Instance()->GetTime() < mEndTime)
    {
        Increment();
        SimulationTime::Instance()->IncrementTimeOneStep();
    }
}

// Explicit instantiation
template class AbstractAngiogenesisSolver<2> ;
template class AbstractAngiogenesisSolver<3> ;
