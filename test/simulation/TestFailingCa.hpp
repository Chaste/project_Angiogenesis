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

#ifndef TESTSPHEROIDWITHANGIOGENESIS_HPP_
#define TESTSPHEROIDWITHANGIOGENESIS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsMesh.hpp"
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "FakePetscSetup.hpp"

class TestSpheroidWithAngiogenesis : public AbstractCellBasedTestSuite
{
public:

    void TestCaBasedSpheroid() throw (Exception)
    {
        // Create a Potts mesh
        unsigned num_x = 31;
        unsigned num_y = 31;
        unsigned num_z = 21;
        PottsMeshGenerator<3> generator(num_x, 0, 0, num_y, 0, 0, num_z, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        double spacing = 40.0; // Fails for this spacing
        //double spacing = 14.0; // Ok for this spacing
        p_mesh->Scale(spacing, spacing, spacing);

        // Create cells in a cylinder
        c_vector<double,3> centre;
        centre[0] = spacing * (num_x - 1)/2.0;
        centre[1] = spacing * (num_y - 1)/2.0;
        centre[2] = 0.0;
        double radius = spacing * 1.2;
        std::vector<unsigned> location_indices;
        for(unsigned kdx=0; kdx<num_z; kdx++)
        {
            for(unsigned jdx=0; jdx<num_y; jdx++)
            {
                for(unsigned idx=0; idx<num_x; idx++)
                {
                    unsigned location_index = idx + num_x * jdx + num_x * num_y * kdx;
                    c_vector<double,3> location;
                    location[0] = double(idx) * spacing;
                    location[1] = double(jdx) * spacing;
                    location[2] = 0.0;

                    double z_location = (double(kdx) * spacing);
                    if(norm_2(location - centre) <= radius && z_location >= 5.0*spacing && z_location <= 25.0 * spacing)
                    {
                        location_indices.push_back(location_index);
                    }
                }
            }
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);

        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestFailingCa");
        simulator.SetDt(1.0);
        simulator.SetEndTime(10.0);
        simulator.Solve();
    }
};

#endif /* TESTSPHEROIDWITHANGIOGENESIS_HPP_ */
