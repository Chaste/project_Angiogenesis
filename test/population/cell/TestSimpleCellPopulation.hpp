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

#ifndef TESTSIMPLECELLPOPULATION_HPP_
#define TESTSIMPLECELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "SmartPointers.hpp"
#include "SimpleCell.hpp"
#include "SimpleCellPopulation.hpp"

class TestSimpleCellPopulation : public CxxTest::TestSuite
{
public:

    void TestGridBasedCellGenerator()
    {
        // Create a helper population
        SimpleCellPopulation population;
        population.GenerateCellsOnGrid(10, 10, 10, 10.0);

        // Write the population to file
        OutputFileHandler output_file_handler("TestSimpleCellPopulation", false);
        population.Write(output_file_handler.GetOutputDirectoryFullPath()+"/cell_population.vtp");
    }

    void TestGridBasedCellGeneratorWithPart()
    {
        // Create the part to be filled with cells
        boost::shared_ptr<Part> p_part = Part::Create();
        boost::shared_ptr<Polygon >p_circle = p_part->AddCircle(50.0, zero_vector<double>(3), 6);
        p_part->Extrude(p_circle, 100.0);

        // Create a helper population
        SimpleCellPopulation population;
        population.GenerateCellsOnGrid(p_part, 10.0);

        // Write the population to file
        OutputFileHandler output_file_handler("TestSimpleCellPopulation", false);
        population.Write(output_file_handler.GetOutputDirectoryFullPath()+"/cell_population_circle.vtp");
        p_part->Write(output_file_handler.GetOutputDirectoryFullPath()+"/circle.vtp");
    }
};

#endif /*TESTSIMPLECELLPOPULATION_HPP_*/
