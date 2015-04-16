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

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartVasculaturePointers.hpp"
#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "VasculatureData.hpp"
#include "ChastePoint.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "FakePetscSetup.hpp"

#include "Debug.hpp"

class TestVascularNode: public AbstractCellBasedTestSuite
{

public:

	void TestConstructor() throw(Exception)
    {
        // Make some nodes and node pointers using explicit coordinates and Chaste points
		VascularNode<2> node1(0.0, 0.0);
        VascularNode<3> node2(ChastePoint<3> (1.0, 1.0, 1.0));
        MAKE_VN_PTR_ARGS(VascularNode<3>, pNode3, (2.0, 2.0, 3.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, pNode4, (ChastePoint<2> (3.0, 4.0)));

        // Test simple Getters and Setters
        node1.SetId(5u);
        std::string label = "Inlet";
        node1.SetLabel(label);
        node1.SetPressure(5.0);
        node1.SetRadius(10.0);

        TS_ASSERT_EQUALS(node1.GetId(), 5u);
        TS_ASSERT_EQUALS(node1.rGetLabel().c_str(), label.c_str());
        TS_ASSERT_DELTA(node1.GetPressure(), 5.0, 1.e-6);
        TS_ASSERT_DELTA(node1.GetRadius(), 10.0, 1.e-6);

        // Test setting location and coincident methods
        ChastePoint<2> point(3.0, 4.0);
        node1.SetLocation(point);
        TS_ASSERT(node1.IsCoincident(point));
        TS_ASSERT(node1.IsCoincident(pNode4));

        // Test distance calculation method
        ChastePoint<2> test_point(2.0, 2.0);
        ChastePoint<3> test_point2(2.0, 2.0, 2.0);
        TS_ASSERT_DELTA(node1.GetDistance(test_point), std::sqrt(5.0), 1.e-6);
        TS_ASSERT_DELTA(node2.GetDistance(test_point2), std::sqrt(3.0), 1.e-6);
    }

	void TestAccessingData() throw(Exception)
    {
        // Make a node
        VascularNode<3> node(1.0, 1.0, 2.0);

        // Set some data
        double radius = 5.5;
        std::string key ="radius";
        node.SetData(key, radius);

        // Check the key is set
        TS_ASSERT(node.HasDataKey("radius"));
        TS_ASSERT_EQUALS(node.GetDataKeys()[0].c_str(), key.c_str());

        bool value_is_castable_to_double = true;
        TS_ASSERT_EQUALS(node.GetDataKeys(value_is_castable_to_double)[0].c_str(), key.c_str());

        // Check the key value is retrieved
        TS_ASSERT_DELTA(node.GetData<double>("radius"), radius, 1.e-6);
        TS_ASSERT_DELTA(node.rGetDataContainer().GetData<double>("radius"), radius, 1.e-6);

        // Replace the existing data container with a new one
        VasculatureData data_container;
        double haematocrit = 7.5;
        data_container.SetData("haematocrit", haematocrit);
        node.SetDataContainer(data_container);
        TS_ASSERT_DELTA(node.GetData<double>("haematocrit"), haematocrit, 1.e-6);
    }

	void TestAddingAndRemovingVesselSegments() throw(Exception)
    {
        // Make some nodes
        MAKE_VN_PTR_ARGS(VascularNode<2>, pNode, (4.0, 3.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, pNode2, (4.0, 5.0));
        MAKE_VN_PTR_ARGS(VascularNode<2>, pNode3, (0.5, 0.6));

        // Make some vessel segments
        MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pVesselSegment, (pNode, pNode2));
        MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pVesselSegment2, (pNode2, pNode3));

        TS_ASSERT(pNode->IsAttachedTo(pVesselSegment));
        TS_ASSERT(!pNode3->IsAttachedTo(pVesselSegment));

        TS_ASSERT_DELTA(pVesselSegment->GetLength(), 2.0, 1.e-6);

        // Check that the vessel segments have been suitably added to the nodes.
        TS_ASSERT_EQUALS(pNode->GetNumberOfSegments(), 1u);
        TS_ASSERT_EQUALS(pNode2->GetNumberOfSegments(), 2u);

        // Check that the segments are correctly retrieved from the node.
        TS_ASSERT(pNode2->IsCoincident(pNode2->GetVesselSegment(0)->GetNode(1)));
        TS_ASSERT(pNode2->IsCoincident(pNode2->GetVesselSegments()[0]->GetNode(1)));
        TS_ASSERT_THROWS_THIS(pNode2->GetVesselSegment(3), "Attempted to access a segment with an out of range index.");

        // Check that the vessel segment connectivity is updated when a node is replaced.
        pVesselSegment2->ReplaceNode(1, pNode);
        TS_ASSERT_EQUALS(pNode->GetNumberOfSegments(), 2u);
        TS_ASSERT_EQUALS(pNode3->GetNumberOfSegments(), 0u);

        // Check that a node can't be replaced with one that's already there
        TS_ASSERT_THROWS_THIS(pVesselSegment2->ReplaceNode(0, pNode), "This segment is already attached to this node.");
    }

    void TestAddingAndRemovingCells() throw(Exception)
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
        TS_ASSERT_THROWS_THIS(node.GetCell(), "A Cell has been requested but none have been assigned to this Node.") ;

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
