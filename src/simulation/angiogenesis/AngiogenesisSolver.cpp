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
#include "UblasIncludes.hpp"
#include "RandomNumberGenerator.hpp"
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "VascularNode.hpp"
#include "AngiogenesisSolver.hpp"

#include "Debug.hpp"

template<unsigned DIM>
AngiogenesisSolver<DIM>::AngiogenesisSolver() :
        mpNetwork(),
        mNodeAnastamosisRadius(5.0),
        mpMigrationRule(),
        mpSproutingRule(),
        mpBoundingDomain(),
        mpFileHandler(),
        mpVesselGrid(),
        mpCellPopulation()
{

}

template<unsigned DIM>
AngiogenesisSolver<DIM>::~AngiogenesisSolver()
{

}

template<unsigned DIM>
boost::shared_ptr<AngiogenesisSolver<DIM> > AngiogenesisSolver<DIM>::Create()
{
    MAKE_PTR(AngiogenesisSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
bool AngiogenesisSolver<DIM>::IsSproutingRuleSet()
{
    return bool(mpSproutingRule);
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetAnastamosisRadius(double radius)
{
    mNodeAnastamosisRadius = radius;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetBoundingDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpBoundingDomain = pDomain;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetCellPopulation(boost::shared_ptr<CaBasedCellPopulationWithVessels<DIM> > cell_population)
{
    mpCellPopulation = cell_population;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetMigrationRule(boost::shared_ptr<AbstractMigrationRule<DIM> > pMigrationRule)
{
    mpMigrationRule = pMigrationRule;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetOutputFileHandler(boost::shared_ptr<OutputFileHandler> pHandler)
{
    mpFileHandler = pHandler;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetSproutingRule(boost::shared_ptr<AbstractSproutingRule<DIM> > pSproutingRule)
{
    mpSproutingRule = pSproutingRule;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetVesselGrid(boost::shared_ptr<RegularGrid<DIM> >pVesselGrid)
{
    mpVesselGrid = pVesselGrid;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetVesselNetwork(boost::shared_ptr<VascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::DoSprouting()
{

    // Get the candidate sprouts and set them as migrating
    std::vector<boost::shared_ptr<VascularNode<DIM> > > candidate_sprouts = mpSproutingRule->GetSprouts(mpNetwork->GetNodes());

    for(unsigned idx=0; idx<candidate_sprouts.size(); idx++)
    {
        candidate_sprouts[idx]->SetIsMigrating(true);
    }

    UpdateNodalPositions(true);

}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::UpdateNodalPositions(bool sprouting)
{

    // Move any nodes marked as migrating, either new sprouts or tips
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > tips;
    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        if(sprouting)
        {
            if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 2)
            {
                tips.push_back(nodes[idx]);
            }
        }
        else
        {
            if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 1)
            {
                tips.push_back(nodes[idx]);
            }
        }
    }

    // Do lattice or off lattice movement
    if (mpVesselGrid)
    {
        std::vector<int> indices = mpMigrationRule->GetIndices(tips);

        for(unsigned idx=0; idx<tips.size();idx++)
        {
            if(indices[idx]>=0)
            {
                if(sprouting)
                {
                    mpNetwork->FormSprout(tips[idx]->GetLocation(), ChastePoint<DIM>(mpVesselGrid->GetLocationOf1dIndex(indices[idx])));
                    tips[idx]->SetIsMigrating(false);
                    mpNetwork->UpdateAll();
                }
                else
                {
                    boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(tips[idx]);
                    p_new_node->SetLocation(mpVesselGrid->GetLocationOf1dIndex(indices[idx]));
                    mpNetwork->ExtendVessel(tips[idx]->GetVesselSegment(0)->GetVessel(), tips[idx], p_new_node);
                    tips[idx]->SetIsMigrating(false);
                    p_new_node->SetIsMigrating(true);
                    mpNetwork->UpdateAll();
                }
            }
            else
            {
                if(sprouting && nodes[idx]->GetNumberOfSegments() == 2)
                {
                    tips[idx]->SetIsMigrating(false);
                }
            }
        }
    }
    else
    {

        mpMigrationRule->SetIsSprouting(sprouting);
        std::vector<c_vector<double, DIM> > movement_vectors = mpMigrationRule->GetDirections(tips);

        for(unsigned idx=0; idx<tips.size();idx++)
        {
            if(norm_2(movement_vectors[idx])>0.0)
            {
                bool do_move = true;
                if(mpBoundingDomain)
                {
                    if(!mpBoundingDomain->IsPointInPart(tips[idx]->GetLocationVector() + movement_vectors[idx]))
                    {
                        do_move = false;
                    }
                }

                if(do_move)
                {
                    if(sprouting)
                    {
                        mpNetwork->FormSprout(tips[idx]->GetLocation(), ChastePoint<DIM>(tips[idx]->GetLocationVector() + movement_vectors[idx]));
                        tips[idx]->SetIsMigrating(false);
                        mpNetwork->UpdateAll();
                    }
                    else
                    {
                        boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(tips[idx]);
                        p_new_node->SetLocation(tips[idx]->GetLocationVector() + movement_vectors[idx]);
                        mpNetwork->ExtendVessel(tips[idx]->GetVesselSegment(0)->GetVessel(), tips[idx], p_new_node);
                        tips[idx]->SetIsMigrating(false);
                        p_new_node->SetIsMigrating(true);
                    }
                }
            }
            else
            {
                if(sprouting && nodes[idx]->GetNumberOfSegments() == 2)
                {
                    tips[idx]->SetIsMigrating(false);
                }
            }
        }
    }

    mpNetwork->UpdateAll();
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::DoAnastamosis()
{
    mpNetwork->UpdateAll();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        // If this is currently a tip
        if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 1)
        {

            if(mpVesselGrid)
            {
                std::vector<std::vector<boost::shared_ptr<VascularNode<DIM> > > > point_node_map = mpVesselGrid->GetPointNodeMap();
                unsigned grid_index = mpVesselGrid->GetNearestGridIndex(nodes[idx]->GetLocationVector());

                if(point_node_map[grid_index].size()>=2)
                {
                    boost::shared_ptr<VascularNode<DIM> > p_merge_node = VascularNode<DIM>::Create(nodes[idx]);

                    if(point_node_map[grid_index][0] == nodes[idx])
                    {
                        p_merge_node = mpNetwork->DivideVessel(point_node_map[grid_index][1]->GetVesselSegment(0)->GetVessel(), nodes[idx]->GetLocation());
                    }
                    else
                    {
                        p_merge_node = mpNetwork->DivideVessel(point_node_map[grid_index][0]->GetVesselSegment(0)->GetVessel(), nodes[idx]->GetLocation());
                    }

                    // Replace the tip node with the merge node
                    p_merge_node->SetIsMigrating(false);
                    if(nodes[idx]->GetVesselSegment(0)->GetNode(0) == nodes[idx])
                    {
                        nodes[idx]->GetVesselSegment(0)->ReplaceNode(0, p_merge_node);
                    }
                    else
                    {
                        nodes[idx]->GetVesselSegment(0)->ReplaceNode(1, p_merge_node);
                    }

                    mpNetwork->UpdateAll();
                }

            }
            else
            {
                // Get the nearest segment and check if it is close enough to the node for a merge
                std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> segment_pair = mpNetwork->GetNearestSegment(nodes[idx], false);

                if(segment_pair.second <= mNodeAnastamosisRadius && nodes[idx]->GetVesselSegment(0)->GetLength() > segment_pair.second)
                {
                    // If there is a non-zero anastamosis radius move the tip onto the segment
                    c_vector<double, DIM> original_location = nodes[idx]->GetLocationVector();
                    if(mNodeAnastamosisRadius > 0.0)
                    {
                        c_vector<double, DIM> divide_location = segment_pair.first->GetPointProjection(original_location, true);
                        nodes[idx]->SetLocation(divide_location);
                    }
                    boost::shared_ptr<VascularNode<DIM> > p_merge_node = mpNetwork->DivideVessel(segment_pair.first->GetVessel(),
                                                                                                 nodes[idx]->GetLocation());

                    // Replace the node at the end of the migrating tip with the merge node
                    if((nodes[idx]->GetVesselSegment(0)->GetNode(0) == p_merge_node) ||
                            (nodes[idx]->GetVesselSegment(0)->GetNode(1) == p_merge_node))
                    {
                        nodes[idx]->SetLocation(original_location);
                    }
                    else
                    {
                        p_merge_node->SetIsMigrating(false);
                        if(nodes[idx]->GetVesselSegment(0)->GetNode(0) == nodes[idx])
                        {
                            nodes[idx]->GetVesselSegment(0)->ReplaceNode(0, p_merge_node);
                        }
                        else
                        {
                            nodes[idx]->GetVesselSegment(0)->ReplaceNode(1, p_merge_node);
                        }
                    }
                    mpNetwork->UpdateAll();
                }
            }
        }
    }
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::Increment()
{
    if(!mpNetwork)
    {
        EXCEPTION("The angiogenesis solver needs an initial vessel network");
    }

    if(mpCellPopulation)
    {
            mpCellPopulation->UpdateVascularCellPopulation();
    }
    else
    {
        // If doing lattice based, add the network to the lattice
        if(mpVesselGrid)
        {
            mpVesselGrid->SetVesselNetwork(mpNetwork);
        }

        // Move any migrating nodes
        UpdateNodalPositions();

        // Check for anastamosis
        DoAnastamosis();

        // Do sprouting
        if(mpSproutingRule)
        {
            mpSproutingRule->SetVesselNetwork(mpNetwork);

            if(mpVesselGrid)
            {
                mpSproutingRule->SetGrid(mpVesselGrid);
            }

            DoSprouting();
            DoAnastamosis();
        }
    }
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::Run(bool writeOutput)
{
    // Loop for the duration of the simulation time
    while(!SimulationTime::Instance()->IsFinished())
    {
        // Write the vessel network if appropriate
        if(writeOutput && mpFileHandler && mpNetwork)
        {
            mpNetwork->Write(mpFileHandler->GetOutputDirectoryFullPath() + "/vessel_network_" +
                             boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()) + ".vtp");
        }

        // Increment the solver and simulation time
        Increment();
        SimulationTime::Instance()->IncrementTimeOneStep();
    }
}

// Explicit instantiation
template class AngiogenesisSolver<2> ;
template class AngiogenesisSolver<3> ;
