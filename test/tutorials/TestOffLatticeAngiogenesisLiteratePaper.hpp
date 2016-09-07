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

#ifndef TESOFFLATTICEANGIOGENESISLITERATEPAPER_HPP_
#define TESOFFLATTICEANGIOGENESISLITERATEPAPER_HPP_

/* = An Off Lattice Angiogenesis Tutorial =
 * This tutorial demonstrates most features of the Angiogenesis Project code. It can be run through
 * to get a rough idea of how the code works, and then individual components can
 * be looked at in more detailed, dedicated tutorials.
 *
 * The following is covered:
 * * Defining domains using geometry primitives
 * * Setting up vessel networks and cell populations
 * * Setting up PDEs and flow problems
 * * Setting up an angiogenesis solver
 * * Interacting with Cell Based Chaste
 *
 * = The Test =
 * Start by introducing the necessary header files. The first contain functionality for setting up unit tests,
 * smart pointer tools and output management.
 */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "RandomNumberGenerator.hpp"
/*
 * Dimensional analysis.
 */
#include "DimensionalChastePoint.hpp"
#include "UnitCollection.hpp"
#include "Owen11Parameters.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
/*
 * Geometry tools.
 */
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
/*
 * Vessel networks.
 */
#include "VesselNode.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
/*
 * Flow.
 */
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
#include "AlarconHaematocritSolver.hpp"
#include "StructuralAdaptationSolver.hpp"
/*
 * Grids and PDEs.
 */
#include "DiscreteContinuumMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "FiniteElementSolver.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
/*
 * Angiogenesis
 */
#include "OffLatticeSproutingRule.hpp"
#include "OffLatticeMigrationRule.hpp"
#include "AngiogenesisSolver.hpp"
/*
 * This should appear last.
 */
#include "PetscSetupAndFinalize.hpp"
class TestOffLatticeAngiogenesisLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * = Test 1 - Angiogenesis With No Cells =
     * In the first example angiogenesis is simulated without interaction with a cell population.
     */
    void TestVesselsOnly() throw(Exception)
    {
        /*
         * We will work in microns.
         */
        units::quantity<unit::length> reference_length(1.0 * unit::microns);
        /*
         * Create  a cylindrical domain, outside of which no vessels can move. Write it to
         * file for visualization.
         */
        units::quantity<unit::length> domain_radius(0.005 * unit::metres);
        units::quantity<unit::length> domain_height(0.001 * unit::metres);
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->SetReferenceLengthScale(reference_length);
        p_domain->AddCylinder(domain_radius, domain_height, DimensionalChastePoint<3>(domain_radius/reference_length,
                                                                                      domain_radius/reference_length, 0.0, reference_length));
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestOffLatticeAngiogenesisLiteratePaper/TestVesselsOnly", false));
        p_domain->Write(p_handler->GetOutputDirectoryFullPath()+"domain.vtp");
        /*
         * Set up a vessel network, which will span the base of the domain.
         */
        units::quantity<unit::length> target_width(1.1*0.01*unit::metres);
        units::quantity<unit::length> target_height(1.1*0.01*unit::metres);
        units::quantity<unit::length> vessel_length(300.0*unit::microns);
        VesselNetworkGenerator<3> network_generator;
        network_generator.SetReferenceLengthScale(reference_length);
        boost::shared_ptr<VesselNetwork<3> > p_network = network_generator.GenerateHexagonalNetwork(target_width,
                                                                                                    target_height,
                                                                                                    vessel_length);
        // Remove any vessel with both nodes outside the domain
        std::vector<boost::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();
        double rad_square = (domain_radius/reference_length)*(domain_radius/reference_length);
        for(unsigned idx=0;idx<vessels.size();idx++)
        {
            double loc_1x = vessels[idx]->GetStartNode()->rGetLocation()[0]-domain_radius/reference_length;
            double loc_1y = vessels[idx]->GetStartNode()->rGetLocation()[1]-domain_radius/reference_length;
            double loc_2x = vessels[idx]->GetEndNode()->rGetLocation()[0]-domain_radius/reference_length;
            double loc_2y = vessels[idx]->GetEndNode()->rGetLocation()[1]-domain_radius/reference_length;
            if((loc_1x*loc_1x+loc_1y*loc_1y)>rad_square && (loc_2x*loc_2x+loc_2y*loc_2y)>rad_square)
            {
                p_network->RemoveVessel(vessels[idx], true);
            }
        }
        // Set a selection of nodes as inlets and outlets
        RandomNumberGenerator::Instance()->Reseed(1101001);
        std::vector<boost::shared_ptr<VesselNode<3> > > nodes = p_network->GetNodes();
        for(unsigned idx=0;idx<nodes.size();idx++)
        {
            if(nodes[idx]->GetNumberOfSegments()==1)
            {
                unsigned rand_num = RandomNumberGenerator::Instance()->randMod(10);
                if(rand_num==0)
                {
                    nodes[idx]->GetFlowProperties()->SetIsInputNode(true);
                    nodes[idx]->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue());
                }
                else if(rand_num<=1)
                {
                    nodes[idx]->GetFlowProperties()->SetIsOutputNode(true);
                    nodes[idx]->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue());
                }
            }
        }
        /*
         * Use the structural adaptation solver to iterate until the flow reaches steady state
         */
        units::quantity<unit::length> vessel_radius(40.0*unit::microns);
        p_network->SetSegmentRadii(vessel_radius);
        units::quantity<unit::dynamic_viscosity> viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
        std::vector<boost::shared_ptr<VesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(viscosity);
        }
        p_network->Write(p_handler->GetOutputDirectoryFullPath()+"initial_network.vtp");

        boost::shared_ptr<DimensionalSimulationTime> p_simulation_time = DimensionalSimulationTime::Instance();
////        SimulationTime::Instance()->SetStartTime(0.0);
//        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<VesselImpedanceCalculator<3> > p_impedance_calculator = boost::shared_ptr<VesselImpedanceCalculator<3> >(new VesselImpedanceCalculator<3>);
//        StructuralAdaptationSolver<3> structural_adaptation_solver;
//        structural_adaptation_solver.SetVesselNetwork(p_network);
//        structural_adaptation_solver.SetWriteOutput(false);
//        structural_adaptation_solver.SetTolerance(0.0001);
//        structural_adaptation_solver.SetMaxIterations(1);
//        structural_adaptation_solver.AddPreFlowSolveCalculator(p_impedance_calculator);
//        structural_adaptation_solver.SetTimeIncrement(0.0001 * unit::seconds);
////        structural_adaptation_solver.Solve();
//
//        p_network->Write(p_handler->GetOutputDirectoryFullPath()+"network_initial_sa.vtp");
        p_simulation_time->Destroy();
    }
};

#endif /*TESTMICROPOCKET_HPP_*/
