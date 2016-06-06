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

#ifndef TESTBUILDVESSELNETWORKLITERATEPAPER_HPP_
#define TESTBUILDVESSELNETWORKLITERATEPAPER_HPP_

/* = Introduction =
 * This tutorial is designed to introduce the C++ interface for modelling vessel networks. An equivalent Python tutorial
 * is [wiki:PaperTutorials/Angiogenesis/PythonBuildVesselNetwork here]. It is advised that you at least read [wiki:UserTutorials/WritingTests
 * the Chaste tutorial on writing tests] before proceeding with this one.
 *
 * This tutorial covers:
 * * Building a network from a collection of nodes, segments and vessels
 * * Writing networks to file and visualizing with Paraview
 * * Building a network using a network generator
 * * Reading a network from file
 *
 * Further functionality is gradually introduced over the course of subsequent tutorials.
 *
 * = The Test =
 * We start by introducing the necessary header files. The first contain functionality for setting up unit tests
 */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
/*
 * Boost shared pointers are used extensively in this component. This header contains some useful
 * pointer MACROs.
 */
#include "SmartPointers.hpp"
/*
 * These headers contain the building-blocks of the vessel networks; nodes, segments, vessels and the network itself
 */
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VascularNetwork.hpp"
/*
 * Tools for generating vessel networks
 */
#include "VasculatureGenerator.hpp"
/*
 * We need to include this when running in serial
 */
#include "FakePetscSetup.hpp"
/*
 * Tutorials are developed as a series of unit tests using the `CxxTest` framework. We make a single test class, which inherits from
 * `AbstractCellBasedWithTimingsTestSuite`. `AbstractCellBasedWithTimingsTestSuite` adds some useful functionality to the default
 * `CxxTest::TestSuite` class, including setting up timers and initializing random number generators.
 */
class TestBuildVesselNetworkLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * = Test 1 - Building a vessel network manually, writing it to file and visualizing it =
     *
     * In the first test we will build a vessel network from its constituent components; nodes, segments and vessels. We will do some
     * simple tests to make sure the network has been formed as expected. Then we write the network to file and visualize it in Paraview.
     */

    void TestBuildNetworkManually() throw (Exception)
    {
        /*
         * First we make some nodes, which are point features from which vessels can be constructed. They are initialized with a location.
         * All vessel network components are created using special factory methods which return shared pointers, rather than being created
         * directly through their constructors. Vessel network components are templated over spatial dimension, and can be 2D or 3D. We will
         * create a Y shaped network. Later we will learn how to build up networks in a more efficient manner.
         */
        double length = 100.0;
        boost::shared_ptr<VascularNode<2> > p_node_1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node_2 = VascularNode<2>::Create(length, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node_3 = VascularNode<2>::Create(2.0 * length, length);
        boost::shared_ptr<VascularNode<2> > p_node_4 = VascularNode<2>::Create(2.0*length, -length);

        /*
         * Next we make vessel segments and vessels. Vessel segments are straight-line features which contain a vascular node at each end. Vessels
         * can be constructed from multiple vessel segments, but in this case each vessel just has a single segment.
         */
        boost::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
        boost::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
        boost::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
        boost::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
        boost::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
        boost::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);

        /*
         * Now we can add our vessels to a vessel network.
         */
        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessel(p_vessel_1);
        p_network->AddVessel(p_vessel_2);
        p_network->AddVessel(p_vessel_3);

        /*
         * We use our test framework to make sure that the network has been created correctly by checking the number of vessels and nodes
         */
        TS_ASSERT_EQUALS(p_network->GetNumberOfNodes(), 4u);
        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(), 3u);

        /*
         * Next we write out network to file. We use the Chaste `OutputFileHandler` functionality to management the output location
         * and the pointer MACRO `MAKE_PTR_ARGS` to quickly make a smart pointer. Networks are written using VTK's PolyData format,
         * which should have a .vtp extension.
         */
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBuildVesselNetworkLiteratePaper"));
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "bifurcating_network.vtp");

        /*
         * Now we can visualize then network in Paraview. See the tutorial [wiki:UserTutorials/VisualizingWithParaview here], to get started. To view the network import the file
         * `TestBuildVesselNetworkLiteratePaper\bifurcating_network.vtp` into Paraview. For a nicer rendering you can do `Filters->Alphabetical->Tube`.
         */
    }

    /*
     * = Test 2 - Building a vessel network using a generator and reading from file =
     *
     * In this test we use a built in generator to automatically construct a network. We then write it to file, read it back in and check
     * that it is restored as expected.
     */
    void TestBuildNetworkFromGeneratorAndReadFromFile() throw (Exception)
    {
        /*
         * We create a hexagonal network in 3D space using a generator. We specify the target network width and height and the desired vessel
         * length.
         */
        VasculatureGenerator<3> network_generator = VasculatureGenerator<3>();
        double target_width = 600.0;
        double target_height = 800.0;
        double length = 100.0;
        boost::shared_ptr<VascularNetwork<3> > p_network = network_generator.GenerateHexagonalNetwork(target_width, target_height, length);

        /*
         * Get the number of nodes and vessels for testing later, and write the network to file as before.
         */
        unsigned number_of_nodes = p_network->GetNumberOfNodes();
        unsigned number_of_vessels = p_network->GetNumberOfVessels();
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestBuildVesselNetworkLiteratePaper"));
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "hexagonal_network.vtp");

        /*
         * We use our generator to read the network back in from the VTK file.
         */
        boost::shared_ptr<VascularNetwork<3> > p_network_from_file =
                network_generator.GenerateNetworkFromVtkFile(p_handler->GetOutputDirectoryFullPath() + "hexagonal_network.vtp");

        /*
         * Finally we check that the network has been correctly read back in using our unit test framework
         */
        TS_ASSERT_EQUALS(p_network_from_file->GetNumberOfNodes(), number_of_nodes);
        TS_ASSERT_EQUALS(p_network_from_file->GetNumberOfVessels(), number_of_vessels);
    }
};

#endif /*TESTBUILDVESSELNETWORKLITERATEPAPER_HPP_*/
