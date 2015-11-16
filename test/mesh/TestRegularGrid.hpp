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

#ifndef TESTREGULARGRID_HPP_
#define TESTREGULARGRID_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RegularGrid.hpp"
#include "RandomNumberGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

class TestRegularGrid : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void dontTestPointPointMapGeneration()
    {
        // Set up a grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        std::vector<unsigned> extents(3);
        extents[0] = 100;
        extents[1] = 100;
        extents[2] = 100;
        p_grid->SetExtents(extents);

        // Set up points
        RandomNumberGenerator::Instance()->Reseed(1000);
        std::vector<c_vector<double, 3> > points(10000);
        for (unsigned idx = 0; idx < 10000; idx++)
        {
            c_vector<double, 3> location;
            location[0] = RandomNumberGenerator::Instance()->ranf() * 100.0;
            location[1] = RandomNumberGenerator::Instance()->ranf() * 100.0;
            location[2] = RandomNumberGenerator::Instance()->ranf() * 100.0;
            points[idx] = location;
        }

        // Get a point-point map
        std::vector<std::vector<unsigned> > map = p_grid->GetPointPointMap(points);

        // Make sure all the points are accounted for
        unsigned sum = 0;
        for (unsigned idx = 0; idx < map.size(); idx++)
        {
            sum += map[idx].size();
        }
        TS_ASSERT_EQUALS(sum, 10000);
    }

    void dontTestPointCellMapGeneration()
    {
        // Set up a grid
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        std::vector<unsigned> extents(3);
        extents[0] = 1000;
        extents[1] = 1000;
        extents[2] = 1;
        p_grid->SetExtents(extents);
        p_grid->SetSpacing(0.333);

        // Set up cells
        HoneycombMeshGenerator generator(100, 100);    // Parameters are: cells across, cells up
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Get a point-cell map
        p_grid->SetCellPopulation(cell_population);
        std::vector<std::vector<CellPtr> > map = p_grid->GetPointCellMap();

        // Make sure all the cells are accounted for
        unsigned sum = 0;
        for (unsigned idx = 0; idx < map.size(); idx++)
        {
            sum += map[idx].size();
        }
        TS_ASSERT_EQUALS(sum, 10000);
    }

    void TestInterpolateGridValues() throw (Exception)
    {
        // Set up a grid
        RandomNumberGenerator::Instance()->Reseed(1000);
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        std::vector<unsigned> extents(3);
        extents[0] = 100;
        extents[1] = 100;
        extents[2] = 100;
        p_grid->SetExtents(extents);
        double spacing = 0.33;
        p_grid->SetSpacing(spacing);

        // Set up a function increasing quadratically from bottom front left to top back right
        std::vector<double> my_grid_func(extents[0] * extents[1] * extents[2]);
        for (unsigned idx = 0; idx < extents[2]; idx++)
        {
            for (unsigned jdx = 0; jdx < extents[1]; jdx++)
            {
                for (unsigned kdx = 0; kdx < extents[0]; kdx++)
                {
                    double value = spacing * spacing * double(kdx * kdx + jdx * jdx + idx * idx);
                    unsigned grid_index = kdx + jdx * extents[0] + idx * extents[0] * extents[1];
                    my_grid_func[grid_index] = value;
                }
            }
        }

        // Set up some sample points
        std::vector<c_vector<double, 3> > points(100);
        for (unsigned idx = 0; idx < 100; idx++)
        {
            c_vector<double, 3> location;
            location[0] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            location[1] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            location[2] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            points[idx] = location;
        }

        // Get the interpolated values
        std::vector<double> interpolated_vals = p_grid->InterpolateGridValues(points, my_grid_func);

        // Get the max error
        double max_error = 0.0;
        for (unsigned idx = 0; idx < interpolated_vals.size(); idx++)
        {
            double analytical = points[idx][0] * points[idx][0] + points[idx][1] * points[idx][1]
                    + points[idx][2] * points[idx][2];
            double error = std::abs((analytical - interpolated_vals[idx]) / analytical);
            if (error > max_error)
            {
                max_error = error;
            }
        }
        TS_ASSERT(max_error < 0.1);
        std::cout << "Max Error: " << max_error << std::endl;
    }

    void TestInterpolateGridValuesWithVtk() throw (Exception)
    {
        // Set up a grid
        RandomNumberGenerator::Instance()->Reseed(1000);
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        std::vector<unsigned> extents(3);
        extents[0] = 100;
        extents[1] = 100;
        extents[2] = 100;
        p_grid->SetExtents(extents);
        double spacing = 0.33;
        p_grid->SetSpacing(spacing);

        // Set up a function increasing quadratically from bottom front left to top back right
        std::vector<double> my_grid_func(extents[0] * extents[1] * extents[2]);
        for (unsigned idx = 0; idx < extents[2]; idx++)
        {
            for (unsigned jdx = 0; jdx < extents[1]; jdx++)
            {
                for (unsigned kdx = 0; kdx < extents[0]; kdx++)
                {
                    double value = spacing * spacing * double(kdx * kdx + jdx * jdx + idx * idx);
                    unsigned grid_index = kdx + jdx * extents[0] + idx * extents[0] * extents[1];
                    my_grid_func[grid_index] = value;
                }
            }
        }

        // Set up some sample points
        std::vector<c_vector<double, 3> > points(100);
        for (unsigned idx = 0; idx < 100; idx++)
        {
            c_vector<double, 3> location;
            location[0] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            location[1] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            location[2] = RandomNumberGenerator::Instance()->ranf() * 30.0;
            points[idx] = location;
        }

        // Get the interpolated values
        std::vector<double> interpolated_vals = p_grid->InterpolateGridValues(points, my_grid_func, true);

        // Get the max error
        double max_error = 0.0;
        for (unsigned idx = 0; idx < interpolated_vals.size(); idx++)
        {
            double analytical = points[idx][0] * points[idx][0] + points[idx][1] * points[idx][1]
                    + points[idx][2] * points[idx][2];
            double error = std::abs((analytical - interpolated_vals[idx]) / analytical);
            if (error > max_error)
            {
                max_error = error;
            }
            std::cout << "vals: " << analytical << "," << interpolated_vals[idx] << std::endl;

        }
        std::cout << "Max Error: " << max_error << std::endl;
        TS_ASSERT(max_error < 0.1);
    }
};

#endif /*TESTREGULARGRID_HPP_*/
