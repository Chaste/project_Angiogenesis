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
#include "SimpleFlowSolver.hpp"
#include "CaVesselSegment.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkProbeFilter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include "Debug.hpp"

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::AbstractAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork) :
        mpNetwork(pNetwork),
        mGrowthVelocity(10.0),
        mTimeIncrement(1.0),
        mEndTime(10.0),
        mOutputFrequency(),
        mOutputDirectory(),
        mNodeAnastamosisRadius(0.0),
        mPdeSolvers(),
        mSolveFlow(false),
        mSproutingProbability(0.0),
        mTimeStep(1.0)
{

}

template<unsigned DIM>
AbstractAngiogenesisSolver<DIM>::~AbstractAngiogenesisSolver()
{

}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetTimeStep(double timeStep)
{
    mTimeStep = timeStep;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > AbstractAngiogenesisSolver<DIM>::GetPdeSolvers()
{
    return mPdeSolvers;
}

template<unsigned DIM>
c_vector<double, DIM> AbstractAngiogenesisSolver<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection)
{
    c_vector<double, DIM> new_direction = currentDirection;

    return new_direction;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetSproutingProbability(double sproutingProbability)
{
    mSproutingProbability = sproutingProbability;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetSolveFlow(bool solveFlow)
{
    mSolveFlow = solveFlow;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::SetOutputDirectory(const std::string& rDirectory)
{
    mOutputDirectory = rDirectory;
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::DoSprouting()
{
    // Randomly sprout if sufficiently far from an existing vessel end node
    mpNetwork->UpdateNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        double prob = RandomNumberGenerator::Instance()->ranf();
        if(nodes[idx]->GetNumberOfSegments()==2 && prob < mSproutingProbability)
        {
            // Do the sprouting. Random normal direction to the existing segment average normals.
            // Get the cross product of the segment tangents
//            c_vector<double, DIM> cross_product = VectorProduct(nodes[idx]->GetVesselSegments()[0]->GetUnitTangent(),
//                                                                nodes[idx]->GetVesselSegments()[1]->GetUnitTangent());
//            double sum = 0.0;
//            for(unsigned jdx=0; jdx<DIM; jdx++)
//            {
//                sum += cross_product[jdx];
//            }
//            if (sum==0.0)
//            {
//                // parallel segments, chose
//            }
            // The sprout will be in the +- x direction, but not along existing vessel directions
            c_vector<double, DIM> sprout_direction;
            if(RandomNumberGenerator::Instance()->ranf()>=0.5)
            {
                sprout_direction = unit_vector<double>(3,0);
            }
            else
            {
                sprout_direction = -unit_vector<double>(3,0);
            }

            // Ensure it is not along the segment vectors
            bool is_along_segment_1 = std::abs(inner_prod(sprout_direction,nodes[idx]->GetVesselSegments()[0]->GetUnitTangent())/
                    (norm_2(sprout_direction)*norm_2(nodes[idx]->GetVesselSegments()[0]->GetUnitTangent()))) > 1 - 1.e-6;
            bool is_along_segment_2 = std::abs(inner_prod(sprout_direction,nodes[idx]->GetVesselSegments()[1]->GetUnitTangent())/
                    (norm_2(sprout_direction)*norm_2(nodes[idx]->GetVesselSegments()[1]->GetUnitTangent()))) > 1 - 1.e-6;

            if(!is_along_segment_1 && !is_along_segment_2)
            {
                mpNetwork->FormSprout(nodes[idx]->GetLocation(), ChastePoint<DIM>(nodes[idx]->GetLocationVector() + mGrowthVelocity*sprout_direction));
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

            // If there is a PDE get the direction of highest solution gradient for the specified species
            if(mPdeSolvers.size()>0)
            {
                int species_index = -1;
                for(unsigned jdx=0; jdx<mPdeSolvers.size(); jdx++)
                {
                    std::string species_name = mPdeSolvers[jdx]->GetPde()->GetVariableName();
                    if(species_name == speciesLabel)
                    {
                        species_index = jdx;
                    }
                }

                if(species_index>=0)
                {
                    // Make points
                    std::vector<c_vector<double, DIM> > locations;
                    locations.push_back(nodes[idx]->GetLocationVector());
                    locations.push_back(locations[0] + mGrowthVelocity * unit_vector<double>(DIM,0));
                    locations.push_back(locations[0] - mGrowthVelocity * unit_vector<double>(DIM,0));
                    locations.push_back(locations[0] + mGrowthVelocity * unit_vector<double>(DIM,1));
                    locations.push_back(locations[0] - mGrowthVelocity * unit_vector<double>(DIM,1));
                    if(DIM==3)
                    {
                        locations.push_back(locations[0] + mGrowthVelocity * unit_vector<double>(DIM,2));
                        locations.push_back(locations[0] - mGrowthVelocity * unit_vector<double>(DIM,2));
                    }

                    vtkSmartPointer<vtkPolyData> p_polydata = vtkSmartPointer<vtkPolyData>::New();
                    vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();
                    p_points->SetNumberOfPoints(locations.size());
                    for(unsigned idx=0; idx< locations.size(); idx++)
                    {
                        if(DIM==3)
                        {
                            p_points->SetPoint(idx, locations[idx][0], locations[idx][1], locations[idx][2]);
                        }
                        else
                        {
                            p_points->SetPoint(idx, locations[idx][0], locations[idx][1], 0.0);
                        }
                    }
                    p_polydata->SetPoints(p_points);

                    vtkSmartPointer<vtkProbeFilter> p_probe_filter = vtkSmartPointer<vtkProbeFilter>::New();
                    p_probe_filter->SetInput(p_polydata);
                    p_probe_filter->SetSource(mPdeSolvers[species_index]->GetSolution());
                    p_probe_filter->Update();
                    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
                }
            }

            // Create a new segment along the growth vector
            boost::shared_ptr<VascularNode<DIM> >  p_new_node = VascularNode<DIM>::Create(nodes[idx]);
            p_new_node->SetLocation(nodes[idx]->GetLocationVector() + mGrowthVelocity * GetGrowthDirection(direction));

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
    mpNetwork->UpdateVesselIds();
    mpNetwork->Write(mOutputDirectory + "/VesselNetwork_inc_" + boost::lexical_cast<std::string>(counter)+".vtp");
    // If there is a flow problem solve it
    SimpleFlowSolver<DIM> flow_solver;
    PoiseuilleImpedanceCalculator<DIM> impedance_calculator;
    if(mSolveFlow)
    {
        mpNetwork->UpdateVesselNodes();
        impedance_calculator.Calculate(mpNetwork);
        flow_solver.SetUp(mpNetwork);
        flow_solver.Implement(mpNetwork);
    }
    // If there is a PDE solve them
    if(mPdeSolvers.size()>0)
    {
        for(unsigned idx=0; idx<mPdeSolvers.size(); idx++)
        {
            mPdeSolvers[idx]->SetVesselNetwork(mpNetwork);
            mPdeSolvers[idx]->SetWorkingDirectory(mOutputDirectory);
            std::string species_name = mPdeSolvers[idx]->GetPde()->GetVariableName();
            mPdeSolvers[idx]->SetFileName("/" + species_name +"_solution_" + boost::lexical_cast<std::string>(counter)+".vti");

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
            mPdeSolvers[idx]->Solve(true);
        }
    }
    while(current_time < mEndTime)
    {
        current_time += mTimeIncrement;

        // Move any migrating nodes
        UpdateNodalPositions();

        DoAnastamosis();

        if(mSproutingProbability > 0.0)
        {
            DoSprouting();
        }

        // Do anastamosis
        DoAnastamosis();
        mpNetwork->MergeCoincidentNodes();
        if(mSolveFlow)
        {
            mpNetwork->UpdateVesselNodes();
            impedance_calculator.Calculate(mpNetwork);
            flow_solver.SetUp(mpNetwork);
            flow_solver.Implement(mpNetwork);
        }

        counter++;
        if(mOutputFrequency > 0 && counter % mOutputFrequency == 0)
        {
            mpNetwork->UpdateVesselIds();
            mpNetwork->Write(mOutputDirectory + "/VesselNetwork_inc_" + boost::lexical_cast<std::string>(counter)+".vtp");
            if(mPdeSolvers.size()>0)
            {
                for(unsigned idx=0; idx<mPdeSolvers.size(); idx++)
                {
                    mPdeSolvers[idx]->SetVesselNetwork(mpNetwork);
                    std::string species_name = mPdeSolvers[idx]->GetPde()->GetVariableName();
                    mPdeSolvers[idx]->SetFileName("/" + species_name +"_solution_" + boost::lexical_cast<std::string>(counter)+".vti");

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
                    mPdeSolvers[idx]->Solve(true);
                }
            }
        }
        else
        {
            if(mPdeSolvers.size()>0)
            {
                for(unsigned idx=0; idx<mPdeSolvers.size(); idx++)
                {

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

                    mPdeSolvers[idx]->Solve(false);
                }
            }
        }
    }

    // destroy the rng
    RandomNumberGenerator::Destroy();
}

template<unsigned DIM>
void AbstractAngiogenesisSolver<DIM>::AddPdeSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pPdeSolver)
{
    mPdeSolvers.push_back(pPdeSolver);
}

// Explicit instantiation
template class AbstractAngiogenesisSolver<2> ;
template class AbstractAngiogenesisSolver<3> ;
