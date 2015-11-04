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

#ifndef TESTDISTANCETRANSFORM_HPP_
#define TESTDISTANCETRANSFORM_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include "Part.hpp"
#include "DistanceTransform.hpp"
#include "CaVascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"

class TestDistanceTransform : public CxxTest::TestSuite
{

public:

    void Test3dBifurcationNetwork()
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        c_vector<double,3> start_position = zero_vector<double>(3);
        start_position[2] = vessel_length;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateBifurcationUnit(vessel_length, start_position);

        // Set up the tissue domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(4.0 * vessel_length, 2.0 * vessel_length, 2.0 * vessel_length);

        // Set up and run the simulation
        DistanceTransform<3> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetGridFromPart(p_domain, 40.0);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestDistanceTransform/Bifurcation3d", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.Setup();
        solver.Solve(true);
    }
};

#endif /*TESTDISTANCETRANSFORM_HPP_*/
