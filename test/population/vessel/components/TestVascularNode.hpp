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

#ifndef TESTVASCULARNODE_HPP_
#define TESTVASCULARNODE_HPP_

#include "Exception.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ChastePoint.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "SmartVasculaturePointers.hpp"
#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "VasculatureData.hpp"
#include "NodeFlowProperties.hpp"
#include "FakePetscSetup.hpp"

class TestVascularNode : public AbstractCellBasedTestSuite
{

public:

    void TestConstructorAndLocationMethods() throw (Exception)
    {
        // Set up some points and locations
        ChastePoint<3> point1(1.0, 2.0, 3.0);
        c_vector<double, 2> location1;
        location1[0] = 1.0;
        location1[1] = 2.0;
        ChastePoint<2> point2(5.0, 6.0);

        // Regular Constructors
        VascularNode<2> node1(0.0, 0.0);
        VascularNode<3> node2(point1);
        VascularNode<2> node3(location1);
        VascularNode<2> node4(node3);

        // Pointer Factory Constructors
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node_1, (2.0, 3.0, 4.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_2, (point2));
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node_3, (location1));
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_4, (node3));
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_5, (p_node_4));

        // Test the location methods
        TS_ASSERT_DELTA(p_node_1->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(p_node_1->GetLocation()[1], 3.0, 1.e-6);
        TS_ASSERT_DELTA(p_node_1->GetLocation()[2], 4.0, 1.e-6);
        TS_ASSERT_DELTA(p_node_5->GetLocationVector()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(p_node_5->GetLocationVector()[1], 2.0, 1.e-6);
    }

    void TestSimpleGetAndSetMethods() throw (Exception)
    {
        // Make a node
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node, (1.0, 2.0, 3.0));

        // Test simple Getters and Setters
        p_node->SetId(5u);
        std::string label = "Inlet";
        p_node->SetLabel(label);
        p_node->GetFlowProperties()->SetPressure(5.0);
        p_node->SetRadius(10.0);
        p_node->GetFlowProperties()->SetIsInputNode(true);
        p_node->GetFlowProperties()->SetIsOutputNode(true);
        p_node->SetIsMigrating(true);

        TS_ASSERT_EQUALS(p_node->GetId(), 5u);
        TS_ASSERT_EQUALS(p_node->rGetLabel().c_str(), label.c_str());
        TS_ASSERT_DELTA(p_node->GetFlowProperties()->GetPressure(), 5.0, 1.e-6);
        TS_ASSERT_DELTA(p_node->GetRadius(), 10.0, 1.e-6);
        TS_ASSERT(p_node->GetFlowProperties()->IsInputNode());
        TS_ASSERT(p_node->GetFlowProperties()->IsOutputNode());
        TS_ASSERT(p_node->IsMigrating());

        // Test setting node flow properties
        NodeFlowProperties node_flow_properties;
        node_flow_properties.SetPressure(12.0);
        p_node->SetFlowProperties(node_flow_properties);

        // Check the data map for the vtk writer
        std::map<std::string, double> vtk_data = p_node->GetVtkData();
        TS_ASSERT_DELTA(vtk_data["Node Id"], 5.0, 1.e-6);
        TS_ASSERT_DELTA(vtk_data["Node Radius"], 10.0, 1.e-6);
        TS_ASSERT_DELTA(vtk_data["Node Pressure"], 12.0, 1.e-6);
        TS_ASSERT_DELTA(vtk_data["Node Is Input"], 0.0, 1.e-6);
        TS_ASSERT_DELTA(vtk_data["Node Is Output"], 0.0, 1.e-6);
        TS_ASSERT_DELTA(vtk_data["Node Is Migrating"], 1.0, 1.e-6);
    }

    void TestDistanceAndConincidentMethods() throw (Exception)
    {
        // Set up some points nodes
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node_1, (1.0, 2.0, 3.0));
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node_2, (1.0, 2.0, 3.0));
        MAKE_VN_PTR_ARGS(VascularNode<3>, p_node_3, (4.0, 5.0, 6.0));
        c_vector<double, 3> location1;
        location1[0] = 6.0;
        location1[1] = 7.0;
        location1[2] = 8.0;
        ChastePoint<3> point1(1.0, 2.0, 3.0);

        // Coincident methods
        TS_ASSERT(p_node_1->IsCoincident(p_node_2));
        TS_ASSERT(p_node_1->IsCoincident(point1));

        // Distance methods
        TS_ASSERT_DELTA(p_node_1->GetDistance(p_node_3), std::sqrt(27.0), 1.e-6);
        TS_ASSERT_DELTA(p_node_1->GetDistance(location1), std::sqrt(75.0), 1.e-6);
        TS_ASSERT_DELTA(p_node_1->GetDistance(point1), 0.0, 1.e-6);
    }

    void TestAccessingData() throw (Exception)
    {
        VascularNode<3> node;
        std::string key ="My Key";
        double value = 5.5;
        node.SetData(key, value);

        // Check the key is set
        TS_ASSERT(node.HasDataKey(key));
        TS_ASSERT_EQUALS(node.GetDataKeys()[0].c_str(), key.c_str());

        bool castable_to_double = true;
        TS_ASSERT_EQUALS(node.GetDataKeys(castable_to_double)[0].c_str(), key.c_str());

        // Check the key value is retrieved
        TS_ASSERT_DELTA(node.GetData<double>(key), value, 1.e-6);
        TS_ASSERT_DELTA(node.rGetDataContainer().GetData<double>(key), value, 1.e-6);

        // Replace the existing data container with a new one
        VasculatureData data_container;
        double new_value = 7.5;
        data_container.SetData("New Key", new_value);
        node.SetDataContainer(data_container);
        TS_ASSERT_DELTA(node.GetData<double>("New Key"), new_value, 1.e-6);
    }

    void TestAddingAndRemovingVesselSegments() throw (Exception)
    {
        // Make some nodes
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_1, (0.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_2, (1.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, p_node_3, (2.0));

        // Make some vessel segments
        MAKE_VN_PTR_ARGS(CaVesselSegment<2>, p_segment1, (p_node_1, p_node_2));
        MAKE_VN_PTR_ARGS(CaVesselSegment<2>, p_segment2, (p_node_2, p_node_3));

        // Check that the vessel segments have been suitably added to the nodes.
        TS_ASSERT(p_node_1->IsAttachedTo(p_segment1));
        TS_ASSERT(!p_node_3->IsAttachedTo(p_segment1));
        TS_ASSERT_EQUALS(p_node_1->GetNumberOfSegments(), 1u);
        TS_ASSERT_EQUALS(p_node_2->GetNumberOfSegments(), 2u);

        // Check that the segments are correctly retrieved from the node.
        TS_ASSERT(p_node_2->IsCoincident(p_node_2->GetVesselSegment(0)->GetNode(1)));
        TS_ASSERT(p_node_2->IsCoincident(p_node_2->GetVesselSegments()[0]->GetNode(1)));
        TS_ASSERT_THROWS_THIS(p_node_2->GetVesselSegment(3), "Attempted to access a segment with an out of range index.");

        // Check that the vessel segment connectivity is updated when a node is replaced.
        p_segment2->ReplaceNode(1, p_node_1);
        TS_ASSERT_EQUALS(p_node_1->GetNumberOfSegments(), 2u);
        TS_ASSERT_EQUALS(p_node_3->GetNumberOfSegments(), 0u);

        // Check that a node can't be replaced with one that's already there
        TS_ASSERT_THROWS_THIS(p_segment2->ReplaceNode(0, p_node_1), "This segment is already attached to this node.");
    }

    void TestAddingAndRemovingCells() throw (Exception)
    {
        // Make a node
        VascularNode<2> node(4.0, 3.0);

        // Create some cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        // Should not be able to add a cell without a cell population
        TS_ASSERT(!node.HasCell());
        TS_ASSERT_THROWS_THIS(node.SetCell(cells[0]), "Attempted to add a Cell without first adding a CellPopulation.");

        // Check that a suitable exception is thrown if the node doesn't have a cell yet
        TS_ASSERT_THROWS_THIS(node.GetCell(), "A Cell has been requested but none have been assigned to this Node.");

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices;
        location_indices.push_back(1);

        // Create a cell population
        MAKE_PTR_ARGS(CaBasedCellPopulation<2>, p_cell_population, (*p_mesh, cells, location_indices));
        node.SetCellPopulation(p_cell_population);

        // Now try adding a cell
        node.SetCell(p_cell_population->rGetCells().front());
        TS_ASSERT(node.HasCell());

        // Verify that the node has moved to the cell's location, also checking GetCell method.
        ChastePoint<2> point2(1.0, 0.0);
        TS_ASSERT(node.IsCoincident(point2));
        TS_ASSERT(point2.IsSamePoint(p_cell_population->GetLocationOfCellCentre(node.GetCell())));

        // Try adding a cell not in the population
        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, 1);
        TS_ASSERT_THROWS_THIS(node.SetCell(cells2[0]), "Attempted to add a Cell that is not in the assigned CellPopulation.");

        // Check that setting a new population removes any old cells
        std::vector<unsigned> location_indices2;
        location_indices2.push_back(2);

        MAKE_PTR_ARGS(CaBasedCellPopulation<2>, p_cell_population2, (*p_mesh, cells2, location_indices2));
        node.SetCellPopulation(p_cell_population2);
        TS_ASSERT(!node.HasCell());

        // Now try removing a cell
        node.SetCell(p_cell_population2->rGetCells().front());
        node.RemoveCell();
        TS_ASSERT(!node.HasCell());
    }
};

#endif /*TESTVASCULARNODE_HPP_*/
