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

#ifndef TESTBLOODFLOWLITERATEPAPER_HPP_
#define TESTBLOODFLOWLITERATEPAPER_HPP_

/*  = Modelling Blood Flow Tutorial =
 * This tutorial demonstrates functionality for modelling blood flow, structural adaptation and vessel
 * regression in a vessel network.
 *
 * This tutorial covers:
 * * Managing parameter values
 * * Running a minimal Poiseuille flow simulation and looking at results
 * * Adding haematocrit
 * * Adding structural adaptation in response to flow
 * * Adding vessel regression in low flow regions
 *
 * = The Test =
 * Start by introducing the necessary header files, explained in previous tutorials.
 */
#include <vector>
#include <cxxtest/TestSuite.h>
#include "Owen11Parameters.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "UnitCollection.hpp"
/*
 * A container for flow information (boundary conditions, pressure values) for nodes.
 */
#include "NodeFlowProperties.hpp"
/*
 * A collection of useful literature parameter values and a way to dump values to
 * file after use.
 */
#include "Owen11Parameters.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
/*
 * The flow and haematocrit solvers, along with neccessary calculators.
 */
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
#include "AlarconHaematocritSolver.hpp"
#include "StructuralAdaptationSolver.hpp"
/*
 * Keep this last.
 */
#include "PetscSetupAndFinalize.hpp"

class TestBloodFlowLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * = Test 1 - Simulating 1d Flow in a Bifurcating Network =
     * [[Image(source:/chaste/projects/Angiogenesis/test/tutorials/images/bifurcation_network_flow.png, 45%, align=center, border=1)]]
     *
     * In the first test we will simulate blood flow in a simple bifurcating vessel network. Subsequent tests will add detail in the form of
     * more complex networks, structural adaptation and vessel regression.
     */
    void TestSimpleFlowProblem() throw (Exception)
    {
        /*
         * We will work in microns
         */
        units::quantity<unit::length> reference_length(1.0 * unit::microns);
        /*
         * First make the network using a generator. Start with a simple unit.
         */
        units::quantity<unit::length> vessel_length = 100.0 * reference_length;
        DimensionalChastePoint<2> start_point(0.0, 0.0);
        VesselNetworkGenerator<2> network_generator;
        network_generator.SetReferenceLengthScale(reference_length);
        boost::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateBifurcationUnit(vessel_length, start_point);
        /*
         * Next, pattern it to make a larger network
         */
        std::vector<unsigned> num_units_per_direction;
        num_units_per_direction.push_back(2);
        num_units_per_direction.push_back(0);
        network_generator.PatternUnitByTranslation(p_network, num_units_per_direction);
        /*
         * Specify which nodes will be the inlets and outlets of the network for the flow problem. This information, as well
         * as all other info related to the flow problem, is contained in a `NodeFlowProperties` instance. Also, set the inlet and
         * outlet pressures in Pa or mmHg. We can mix and match these thanks to the dimensional analysis functionality. Because
         * the network is simple, we can figure out which node is which just from their index. For more complicated networks
         * we need to use spatial locators and `NearestNode` type methods.
         */
        units::quantity<unit::pressure> inlet_pressure(50.0 * unit::mmHg);
        p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
        p_network->GetNode(0)->GetFlowProperties()->SetPressure(inlet_pressure);
        /*
         * It would be useful if we had a record of which parameter values we have used in the simulation and where they are sourced in the
         * literature. Instead of manually entering parameter values like above we can use a `ParameterCollection` singleton which allows for
         * some extra metadata storage. We will take some parameter values from a paper by Owen et al. (2011). We add the parameter value
         * to our parameter collection with the annotation 'User', to clarify it is not added automatically by some solver.
         */
        p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
        p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue());
        /*
         * Now set the segment radii and viscosity values.
         */
        units::quantity<unit::length> vessel_radius(GenericParameters::mpCapillaryRadius->GetValue());
        p_network->SetSegmentRadii(vessel_radius);

        units::quantity<unit::dynamic_viscosity> viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
        std::vector<boost::shared_ptr<VesselSegment<2> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(viscosity);
        }
        /*
         * We use a calculator to work out the impedance of each vessel based on assumptions of Poiseuille flow and cylindrical vessels. This
         * updates the value of the impedance in the vessel.
         */
        VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
        impedance_calculator.SetVesselNetwork(p_network);
        impedance_calculator.Calculate();
        /*
         * Check that the impedance is as expected in one of the vessels
         */
        units::quantity<unit::flow_impedance> expected_impedance = 8.0 * viscosity* vessel_length/(M_PI*units::pow<4>(vessel_radius));
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetSegment(0)->GetFlowProperties()->GetImpedance().value(), expected_impedance.value(), 1.e-6);
        /*
         * Now we can solve for the flow rates in each vessel based on the inlet and outlet pressures and impedances. The solver
         * updates the value of pressures and flow rates in each vessel and node in the network.
         */
        FlowSolver<2> flow_solver = FlowSolver<2>();
        flow_solver.SetVesselNetwork(p_network);
        flow_solver.Solve();
        /*
         * Check the pressures, it is expected to drop linearly so should be the average of the input and output half way along the network.
         */
        units::quantity<unit::pressure> expected_pressure = (inlet_pressure + Owen11Parameters::mpOutletPressure->GetValue())/2.0;
        TS_ASSERT_DELTA(p_network->GetNode(7)->GetFlowProperties()->GetPressure().value(), expected_pressure.value(), 1.e-6);
        /*
         * Next we write out the network, including updated flow data, to file.
         */
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBloodFlowLiteratePaper/TestSimpleFlowProblem"));
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "bifurcating_network_results.vtp");
        /*
         * Now we can visualize the results in Paraview. See [wiki:UserTutorials/VisualizingWithParaview here] to get started. To view the network import the file
         * `TestBloodFlowLiteratePaper\bifurcating_network.vtp` into Paraview. For a nicer rendering you can do `Filters->Alphabetical->Tube`.
         * Finally, dump our parameter collection to an xml file and, importantly, clear it for the next test.
         */
        ParameterCollection::Instance()->DumpToFile(p_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
        ParameterCollection::Instance()->Destroy();
    }
    /*
     * = Test 2 - Simulating Haematocrit Transport in 3D =
     * [[Image(source:/chaste/projects/Angiogenesis/test/tutorials/images/haematocrit.png, 25%, align=center, border=1)]]
     *
     * In this test we will simulate haematocrit transport in a 3d vessel network.
     */
    void TestFlowProblemWithHaematocrit() throw (Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBloodFlowLiteratePaper/TestFlowProblemWithHaematocrit", false));
        /*
         * This time we solve a flow problem and then use the solution to calculate the haematocrit distribution,
         * assuming it has no effect on the flow. Set up the network as before
         */
        units::quantity<unit::length> cell_width(25.0 * unit::microns);
        units::quantity<unit::length> target_width = 100.0 * cell_width;
        units::quantity<unit::length> target_height = 30.0 * cell_width;
        units::quantity<unit::length> vessel_length = 4.0 * cell_width;
        VesselNetworkGenerator<3> network_generator;
        network_generator.SetReferenceLengthScale(cell_width);
        boost::shared_ptr<VesselNetwork<3> > p_network = network_generator.GenerateHexagonalNetwork(target_width,
                                                                                                    target_height,
                                                                                                    vessel_length);
        /*
        * We will use a locator to mark the bottom left and top right nodes as respective inlets and outlets
        */
        DimensionalChastePoint<3> inlet_locator(0.0, 0.0, 0.0, cell_width);
        DimensionalChastePoint<3> outlet_locator(target_width/cell_width, target_height/cell_width, 0.0, cell_width);
        boost::shared_ptr<VesselNode<3> > p_inlet_node = p_network->GetNearestNode(inlet_locator);
        boost::shared_ptr<VesselNode<3> > p_outlet_node = p_network->GetNearestNode(outlet_locator);
        p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
        p_inlet_node->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue());
        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
        p_outlet_node->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue());

        /*
         * Now set the segment radii and viscosity values.
         */
        units::quantity<unit::length> vessel_radius(GenericParameters::mpCapillaryRadius->GetValue());
        p_network->SetSegmentRadii(vessel_radius);
        units::quantity<unit::dynamic_viscosity> viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
        std::vector<boost::shared_ptr<VesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(viscosity);
        }
        /*
         * Next some simple extra functionality is demonstrated by mapping the network onto a hemisphere
         */
        std::vector<std::pair<double, double> > extents = p_network->GetExtents();
        double sphere_radius = 400.0;
        double sphere_thickess = 3.0;
        double sphere_azimuth = M_PI;
        double sphere_polar = M_PI/2.0;
        std::vector<boost::shared_ptr<VesselNode<3> > > nodes = p_network->GetNodes();
        for(unsigned idx =0; idx<nodes.size(); idx++)
        {
            double x_frac = (nodes[idx]->rGetLocation()[0]) / (extents[0].second - extents[0].first);
            double azimuth_angle = x_frac * sphere_azimuth;

            double y_frac = (nodes[idx]->rGetLocation()[1] + 3.0) / (extents[1].second - extents[1].first);
            double polar_angle = y_frac * sphere_polar;

            double radius = sphere_radius;
            if(extents[2].second - extents[2].first>0.0)
            {
                double z_frac = nodes[idx]->rGetLocation()[2] / (extents[2].second - extents[2].first);
                radius = sphere_radius - sphere_thickess * z_frac;
            }
            DimensionalChastePoint<3>new_position(radius * std::cos(azimuth_angle) * std::sin(polar_angle),
                                                  radius * std::cos(polar_angle),
                                                  radius * std::sin(azimuth_angle) * std::sin(polar_angle),
                                                  cell_width);
            nodes[idx]->SetLocation(new_position);
        }
        /*
         * Get the impedance
         */
        VesselImpedanceCalculator<3> impedance_calculator = VesselImpedanceCalculator<3>();
        impedance_calculator.SetVesselNetwork(p_network);
        impedance_calculator.Calculate();
        /*
         * Solve the flow problem
         */
        FlowSolver<3> flow_solver = FlowSolver<3>();
        flow_solver.SetVesselNetwork(p_network);
        flow_solver.Solve();
        /*
         * Solve the haematocrit problem
         */
        AlarconHaematocritSolver<3> haematocrit_solver = AlarconHaematocritSolver<3>();
        haematocrit_solver.SetVesselNetwork(p_network);
        haematocrit_solver.SetHaematocrit(0.45);
        haematocrit_solver.Calculate();
        /*
         * Next we write out the network, including updated flow data, to file.
         */
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "network_haematocrit.vtp");
        /*
         * Now we can visualize the results in Paraview. See [wiki:UserTutorials/VisualizingWithParaview here] to get started. To view the network import the file
         * `TestBloodFlowLiteratePaper\bifurcating_network.vtp` into Paraview. For a nicer rendering you can do `Filters->Alphabetical->Tube`.
         */
        ParameterCollection::Instance()->DumpToFile(p_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
        ParameterCollection::Instance()->Destroy();
    }
};
#endif /*TESTBLOODFLOWLITERATEPAPER_HPP_*/
