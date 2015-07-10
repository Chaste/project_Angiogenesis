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

#ifndef TESTVORONOIGENERATOR_HPP_
#define TESTVORONOIGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "Part.hpp"
#include "VoronoiGenerator.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "VasculatureGenerator.hpp"


class TestVoronoiGenerator : public CxxTest::TestSuite
{
public:

    void DontTestSquare()
    {
        boost::shared_ptr<Part> p_part = Part::Create();
        p_part->AddCuboid(3000, 1500, 200);
        VoronoiGenerator generator;
        boost::shared_ptr<Part> p_tesselation = generator.Generate(p_part, std::vector<boost::shared_ptr<Vertex> >(), 1600);
        std::vector<boost::shared_ptr<Vertex> > vertices = p_tesselation->GetVertices();
        for(unsigned idx=0; idx < vertices.size(); idx++)
        {
            vertices[idx]->SetCoordinate(1, 2.0 * vertices[idx]->rGetLocation()[1]);
        }
        OutputFileHandler output_file_handler("TestVoronoiNetwork", false);
        p_tesselation->Write(output_file_handler.GetOutputDirectoryFullPath() + "part.vtp");
    }

    void TestSquareWithPde()
    {
        boost::shared_ptr<Part> p_part = Part::Create();
        p_part->AddCuboid(3000, 1500, 200);
        VoronoiGenerator generator;
        boost::shared_ptr<Part> p_tesselation = generator.Generate(p_part, std::vector<boost::shared_ptr<Vertex> >(), 1000);
        std::vector<boost::shared_ptr<Vertex> > vertices = p_tesselation->GetVertices();
        for(unsigned idx=0; idx < vertices.size(); idx++)
        {
            vertices[idx]->SetCoordinate(1, 2.0 * vertices[idx]->rGetLocation()[1]);
        }
        boost::shared_ptr<CaVascularNetwork<3> >p_network =  p_tesselation->GenerateVesselNetwork();

        // Set up the PDE domain
        boost::shared_ptr<Part> p_domain = Part::Create();
        p_domain->AddCuboid(3000, 3000, 200);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the solver
        FiniteDifferenceSolver solver;
        solver.SetVesselNetwork(p_network);
        solver.SetExtents(p_domain, 20.0);
        solver.SetPde(p_pde);

        OutputFileHandler output_file_handler("TestVoronoiNetwork/Pde", false);
        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
        solver.Solve(true);
    }

    void DontTestBioWithPde()
    {
        FileFinder file_finder("projects/JamesG/test/data/bio.vtp", RelativeTo::ChasteSourceRoot);
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> >p_network =  generator.GenerateNetworkFromVtkFile(file_finder.GetAbsolutePath());

        // Set up the PDE domain
        boost::shared_ptr<Part> p_domain = Part::Create();
        p_domain->AddCuboid(3400, 3400, 20);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the solver
        FiniteDifferenceSolver solver;
        solver.SetVesselNetwork(p_network);
        solver.SetExtents(p_domain, 20.0);
        solver.SetPde(p_pde);

        OutputFileHandler output_file_handler("TestVoronoiNetwork/Bio", false);
        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
        solver.Solve(true);
    }

};

#endif /*TESTVORONOIGENERATOR_HPP_*/
