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
#include "AbstractAngiogenesisSolver.hpp"

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::AbstractAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork, const std::string& rOutputDirectory) :
        mpNetwork(pNetwork),
        mGrowthVelocity(10.0),
        mTimeIncrement(1.0),
        mEndTime(10.0),
        mOutputFrequency(1),
        mOutputDirectory(rOutputDirectory),
        mNodeAnastamosisRadius(0.0),
        mpPdeSolver()
{

}

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::~AbstractAngiogenesisSolver()
{

}

template<unsigned DIM>
c_vector<double, DIM> AbstractAngiogenesisSolver<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection)
{
    c_vector<double, DIM> new_direction = currentDirection;

    return new_direction;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::UpdateNodalPositions()
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
         boost::shared_ptr<VascularNode<DIM> >  p_new_node = VascularNode<DIM>::Create(nodes[idx]);
         p_new_node->SetLocation(nodes[idx]->GetLocationVector() + mGrowthVelocity * GetGrowthDirection(direction));

         if(nodes[idx]->GetVesselSegment(0)->GetVessel()->GetStartNode() == nodes[idx])
         {
             nodes[idx]->GetVesselSegment(0)->GetVessel()->AddSegment(CaVesselSegment<DIM>::Create(p_new_node, nodes[idx]));
         }
         else
         {
             nodes[idx]->GetVesselSegment(0)->GetVessel()->AddSegment(CaVesselSegment<DIM>::Create(nodes[idx], p_new_node));
         }
         nodes[idx]->SetIsMigrating(false);
         p_new_node->SetIsMigrating(true);
     }
    }
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
                remaining_nodes[idx]->SetLocation(divide_location);
                remaining_nodes[idx]->SetIsMigrating(false);
            }
        }
    }
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::Run()
{
    // Loop over the time (replace with simulation time)
    double current_time = 0.0;

    // make a new rng
    RandomNumberGenerator::Instance();

    unsigned counter = 0;
    mpNetwork->MergeCoincidentNodes();
    mpNetwork->Write(mOutputDirectory + "/VesselNetwork_inc_" + boost::lexical_cast<std::string>(counter)+".vtp");

    // If there is a PDE solve it
    if(mpPdeSolver)
    {
        mpPdeSolver->SetVesselNetwork(mpNetwork);
        mpPdeSolver->SetWorkingDirectory(mOutputDirectory);
        mpPdeSolver->SetFileName("/solution_" + boost::lexical_cast<std::string>(counter)+".vti");
        mpPdeSolver->Solve(true);
    }

    while(current_time < mEndTime)
    {
        current_time += mTimeIncrement;

        // Move any migrating nodes
        UpdateNodalPositions();

        // Do anastamosis
        DoAnastamosis();

        mpNetwork->MergeCoincidentNodes();
        counter++;
        if(mOutputFrequency > 0 && counter % mOutputFrequency == 0)
        {
            mpNetwork->Write(mOutputDirectory + "/VesselNetwork_inc_" + boost::lexical_cast<std::string>(counter)+".vtp");
            if(mpPdeSolver)
            {
                mpPdeSolver->SetVesselNetwork(mpNetwork);
                mpPdeSolver->SetFileName("/solution_" + boost::lexical_cast<std::string>(counter)+".vti");
                mpPdeSolver->Solve(true);
            }
        }
        else
        {
            if(mpPdeSolver)
            {
                mpPdeSolver->Solve(false);
            }
        }
    }

    // destroy the rng
    RandomNumberGenerator::Destroy();
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetPdeSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pPdeSolver)
{
    mpPdeSolver = pPdeSolver;

}

// Explicit instantiation
template class AbstractAngiogenesisSolver<2> ;
template class AbstractAngiogenesisSolver<3> ;
