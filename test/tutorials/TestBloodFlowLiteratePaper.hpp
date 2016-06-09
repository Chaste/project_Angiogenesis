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
 * regression using the C++ framework.
 *
 * This tutorial covers:
 * * Running a minimal Poiseuille flow simulation and looking at results
 * * Adding haematocrit
 * * Adding structural adaptation in response to flow
 * * Adding vessel regression in low flow regions
 *
 * = The Test =
 * We start by introducing the necessary header files.
 */
#include <vector>
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
/*
 * The flow solver, individual calculators and the structural adaptation solver
 */
#include "FlowSolver.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
/*
 * We need to include this when running in serial
 */
#include "FakePetscSetup.hpp"

class TestBloodFlowLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * = Test 1 - Simulating 1d flow in a bifurcating network=
     *
     * In the first test we will simulate blood flow in a simple bifurcating vessel network. Subsequent tests will add detail in the form of
     * more complex networks, structural adaptation and vessel regression.
     */
    void TestSimpleFlowProblem() throw (Exception)
    {
        /*
         * First we make the network. It has the same  Y shape as the previous tutorial, but is constructed a little more efficiently.
         * Note that segments do not need to be explicitly created if there is only one per vessel.
         */
        double length = 100.0;
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        nodes.push_back(VascularNode<2>::Create(0.0, 0.0));
        nodes.push_back(VascularNode<2>::Create(length, 0.0));
        nodes.push_back(VascularNode<2>::Create(2.0*length, length));
        nodes.push_back(VascularNode<2>::Create(2.0*length, -length));
        std::vector<boost::shared_ptr<Vessel<2> > > vessels;
        vessels.push_back(Vessel<2>::Create(nodes[0], nodes[1]));
        vessels.push_back(Vessel<2>::Create(nodes[1], nodes[2]));
        vessels.push_back(Vessel<2>::Create(nodes[1], nodes[3]));
        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessels(vessels);
        /*
         * We specify which nodes will be the inlets and outlets of the network for the flow problem. This information, as well
         * as all other info related to the flow problem, is contained in a `NodeFlowProperties` instance. Then we set the inlet and
         * outlet pressures in Pa. Finally, we specify the radius and viscosity of each segment.
         * This is used in the calculation of the impedance of the vessel in the flow problem.
         */
        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(5000.0);
        nodes[2]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[2]->GetFlowProperties()->SetPressure(3000.0);
        nodes[3]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[3]->GetFlowProperties()->SetPressure(3000.0);
        p_network->SetSegmentRadii(10.0);
        std::vector<boost::shared_ptr<VesselSegment<2> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }
        /*
         * We use a calculator to work out the impedance of each vessel based on assumptions of Poiseuille flow and cylindrical vessels. This
         * updates the value of the impedance in the vessel.
         */
        PoiseuilleImpedanceCalculator<2> impedance_calculator = PoiseuilleImpedanceCalculator<2>();
        impedance_calculator.Calculate(p_network);
        /*
         * Now we can solve for the flow rates in each vessel based on the inlet and outlet pressures and impedances. The solver
         * updates the value of pressures and flow rates in each vessel and node in the network.
         */
        FlowSolver<2> flow_solver = FlowSolver<2>();
        flow_solver.SetVesselNetwork(p_network);
        flow_solver.Solve();
        /*
         * We can check to see if the final solution is reasonable
         */
//        TS_ASSERT_EQUALS(p_network->GetNumberOfNodes(), 4u);
//        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(), 3u);

        /*
         * Next we write out the network, including updated flow data, to file.
         */
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBloodFlowLiteratePaper"));
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "bifurcating_network_results.vtp");
        /*
         * Now we can visualize the results in Paraview. See [wiki:UserTutorials/VisualizingWithParaview here] to get started. To view the network import the file
         * `TestBloodFlowLiteratePaper\bifurcating_network.vtp` into Paraview. For a nicer rendering you can do `Filters->Alphabetical->Tube`.
         */
    }

    void TestFlowProblemWithHaematocrit() throw (Exception)
    {
        /*
         * This time we solve a flow problem with haematocrit. We set up the problem as before.
         */
        double length = 100.0;
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        nodes.push_back(VascularNode<2>::Create(0.0, 0.0));
        nodes.push_back(VascularNode<2>::Create(length, 0.0));
        nodes.push_back(VascularNode<2>::Create(2.0*length, length));
        nodes.push_back(VascularNode<2>::Create(2.0*length, -length));
        std::vector<boost::shared_ptr<Vessel<2> > > vessels;
        vessels.push_back(Vessel<2>::Create(nodes[0], nodes[1]));
        vessels.push_back(Vessel<2>::Create(nodes[1], nodes[2]));
        vessels.push_back(Vessel<2>::Create(nodes[1], nodes[3]));
        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessels(vessels);
        /*
         * We specify which nodes will be the inlets and outlets of the network for the flow problem. This information, as well
         * as all other info related to the flow problem, is contained in a `NodeFlowProperties` instance. Then we set the inlet and
         * outlet pressures in Pa. Finally, we specify the radius and viscosity of each segment.
         * This is used in the calculation of the impedance of the vessel in the flow problem.
         */
        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(5000.0);
        nodes[2]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[2]->GetFlowProperties()->SetPressure(3000.0);
        nodes[3]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[3]->GetFlowProperties()->SetPressure(3000.0);
        p_network->SetSegmentRadii(10.0);
        std::vector<boost::shared_ptr<VesselSegment<2> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }
        /*
         * We use a calculator to work out the impedance of each vessel based on assumptions of Poiseuille flow and cylindrical vessels. This
         * updates the value of the impedance in the vessel.
         */
        PoiseuilleImpedanceCalculator<2> impedance_calculator = PoiseuilleImpedanceCalculator<2>();
        impedance_calculator.Calculate(p_network);
        /*
         * Now we can solve for the flow rates in each vessel based on the inlet and outlet pressures and impedances. The solver
         * updates the value of pressures and flow rates in each vessel and node in the network.
         */
        FlowSolver<2> flow_solver = FlowSolver<2>();
        flow_solver.SetVesselNetwork(p_network);
        flow_solver.Solve();
        /*
         * We can check to see if the final solution is reasonable
         */
//        TS_ASSERT_EQUALS(p_network->GetNumberOfNodes(), 4u);
//        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(), 3u);

        /*
         * Next we write out the network, including updated flow data, to file.
         */
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBloodFlowLiteratePaper"));
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "bifurcating_network_results.vtp");
        /*
         * Now we can visualize the results in Paraview. See [wiki:UserTutorials/VisualizingWithParaview here] to get started. To view the network import the file
         * `TestBloodFlowLiteratePaper\bifurcating_network.vtp` into Paraview. For a nicer rendering you can do `Filters->Alphabetical->Tube`.
         */
    }


};

#endif /*TESTBLOODFLOWLITERATEPAPER_HPP_*/
