/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TestOnLatticeRwGrowthDirectionModifier_hpp
#define TestOnLatticeRwGrowthDirectionModifier_hpp

#include <cxxtest/TestSuite.h>
#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OnLatticeRwGrowthDirectionModifier.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "AngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestOnLatticeRwGrowthDirectionModifier : public AbstractCellBasedTestSuite
{

public:

    void TestSingleVesselGrowth() throw(Exception)
    {
        // Make a network
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(100.0, 0.0, 0.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        // Set the one of the nodes to migrate
        p_network->GetVessel(0)->GetEndNode()->SetIsMigrating(true);

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/SingleVesselGrowth/");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiVessel() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 0.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 100.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->FormSprout(ChastePoint<3>(20.0, 0.0, 0.0), ChastePoint<3>(20.0, 10.0, 0.0));
        p_network->FormSprout(ChastePoint<3>(70.0, 100.0, 0.0), ChastePoint<3>(70.0, 90.0, 0.0));

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/Multisegment/");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSprout() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 0.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=1; idx<9; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 0.0, 0.0), ChastePoint<3>(double(idx)*10, 10.0, 0.0));
        }

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/MultiSprout/");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSproutWithPde() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 90.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        for(unsigned idx=1; idx<10; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        // Set up the PDE domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(100, 100, 10);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up the boundary condition
        boost::shared_ptr<DirichletBoundaryCondition<3> > p_vessel_ox_boundary_condition = DirichletBoundaryCondition<3>::Create();
        p_vessel_ox_boundary_condition->SetValue(40.0);
        p_vessel_ox_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
        p_vessel_ox_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);
        p_vessel_ox_boundary_condition->SetVesselNetwork(p_network);

        // Set up and run the solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_solver = FiniteDifferenceSolver<3>::Create();
        p_solver->SetGridFromPart(p_domain, 10.0);
        p_solver->SetPde(p_pde);

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/MultiSproutPde/");
        angiogenesis_solver.AddPdeSolver(p_solver);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSproutWithFlow() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[10]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[10]->GetFlowProperties()->SetPressure(1000);

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 90.0, 0.0));
        }
        top_nodes[10]->GetFlowProperties()->SetIsInputNode(true);
        top_nodes[10]->GetFlowProperties()->SetPressure(3000);
        top_nodes[0]->GetFlowProperties()->SetIsOutputNode(true);
        top_nodes[0]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<10; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        p_network->UpdateSegments();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/MultiSproutFlow/");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestSproutingWithFlow() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<9; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<6; idx+=2)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        p_network->UpdateSegments();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/SproutingFlow/");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }
};

#endif
