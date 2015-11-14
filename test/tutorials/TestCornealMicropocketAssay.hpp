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
 * Redistributions in binary form must reproduce the abovea copyright notice,
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

#ifndef TESTMICROPOCKET_HPP_
#define TESTMICROPOCKET_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include <math.h>

#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OffLatticePrwGrowthDirectionModifier.hpp"
#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OffLatticeSolutionDependentGrowthDirectionModifier.hpp"
#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OffLatticeTipAttractionGrowthDirectionModifier.hpp"
#include "../../src/angiogenesis_solvers/sprouting_rules/OffLatticeRandomNormalSproutingRule.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
#include "PlcMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "FiniteElementSolver.hpp"
#include "DirichletBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AngiogenesisSolver.hpp"
#include "Debug.hpp"

class TestMicropocketGeometry : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestTheGeometry() throw(Exception)
    {
        // Make the domain
        double sphere_radius = 1500.0;
        double sphere_thickness = 100.0;
        double sphere_azimuth_angle = M_PI;
        double sphere_polar_angle = M_PI/2.0;

        MappableGridGenerator generator;
        boost::shared_ptr<Part<3> > p_part = generator.GenerateHemisphere(sphere_radius, sphere_thickness , 10, 10, sphere_azimuth_angle, sphere_polar_angle);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 1.e6);

        // Make the vessel network
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();

        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        double node_spacing = 10.0;
        unsigned num_nodes = 40;
        double frequency = 3.0;
        double amplitude = 20.0;
        double total_length = double(num_nodes-1) * node_spacing;

        for(unsigned idx=0; idx<num_nodes; idx++)
        {
            double x_position = double(idx) * node_spacing;
            double y_position = amplitude * std::sin((x_position / total_length) * 2.0 * M_PI * frequency);
            nodes.push_back(VascularNode<3>::Create(x_position, y_position + 2.0 * amplitude, sphere_thickness / 2.0));
        }

        // Map the nodes to the sphere
        for(unsigned idx =0; idx<nodes.size(); idx++)
        {
            double x_frac = nodes[idx]->GetLocationVector()[0] / total_length;
            double azimuth_angle = x_frac * sphere_azimuth_angle;

            double y_frac = (1.0 - nodes[idx]->GetLocationVector()[1] / (0.5*M_PI*sphere_radius));
            double polar_angle = y_frac * sphere_polar_angle;

            double z_frac = nodes[idx]->GetLocationVector()[2] / (sphere_thickness);
            double radius = sphere_radius - sphere_thickness * z_frac;

            // Get the new x
            c_vector<double, 3> new_position;
            new_position[0] = radius * std::cos(azimuth_angle) * std::sin(polar_angle);

            // Get the new y
            new_position[1] = radius * std::cos(polar_angle);

            // Get the new z
            new_position[2] = radius * std::sin(azimuth_angle) * std::sin(polar_angle);

            nodes[idx]->SetLocation(new_position);
        }

        boost::shared_ptr<CaVessel<3> > p_vessel_1 = CaVessel<3>::Create(nodes);

        // Set up flow properties
        p_network->AddVessel(p_vessel_1);
        p_network->GetVessel(0)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_network->GetVessel(0)->GetStartNode()->GetFlowProperties()->SetPressure(3000.0);
        p_network->GetVessel(0)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_network->GetVessel(0)->GetEndNode()->GetFlowProperties()->SetPressure(1000.0);
        p_network->UpdateSegments();
        p_network->SetSegmentRadii(10.0);
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-9);
        }

        boost::shared_ptr<HybridLinearEllipticPde<3> > p_vegf_pde = HybridLinearEllipticPde<3>::Create();
        p_vegf_pde->SetDiffusionConstant(0.0033);
        p_vegf_pde->SetVariableName("vegf");
        p_vegf_pde->SetLinearInUTerm(-1.e-8);

        boost::shared_ptr<FiniteElementSolver<3> > p_vegf_solver = FiniteElementSolver<3>::Create();

        // Generate the boundary condition region
        boost::shared_ptr<Part<3> > p_boundary_part = Part<3>::Create();
        p_boundary_part->AddCuboid(200.0, 200.0, 3000.0);
        c_vector<double, 3> translation_vector;
        translation_vector[0] = -100.0;
        translation_vector[1] = 1100.0;
        translation_vector[2] = 0.0;
        p_boundary_part->Translate(translation_vector);
        boost::shared_ptr<DirichletBoundaryCondition<3> > p_vegf_boundary_condition = DirichletBoundaryCondition<3>::Create();
        p_vegf_boundary_condition->SetValue(40.0);
        p_vegf_boundary_condition->SetType(BoundaryConditionType::IN_PART);
        p_vegf_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);
        p_vegf_boundary_condition->SetDomain(p_boundary_part);

        p_vegf_solver->AddDirichletBoundaryCondition(p_vegf_boundary_condition);
        p_vegf_solver->SetMesh(p_mesh);
        p_vegf_solver->SetPde(p_vegf_pde);

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);

        boost::shared_ptr<OffLatticeSolutionDependentGrowthDirectionModifier<3> > p_grow_direction_modifier3 = OffLatticeSolutionDependentGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier3->SetSolver(p_vegf_solver);
        p_grow_direction_modifier3->SetStrength(0.5);

        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.05);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestMicropocket");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier2);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier3);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetAnastamosisRadius(5.0);
        angiogenesis_solver.AddPdeSolver(p_vegf_solver);
        angiogenesis_solver.SetBoundingDomain(p_part);
        angiogenesis_solver.SetEndTime(10);
        angiogenesis_solver.Run();
    }
};

#endif /*TESTMICROPOCKET_HPP_*/
