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

#ifndef TESTMAPPABLEGRIDGENERATOR_HPP_
#define TESTMAPPABLEGRIDGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include <math.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
#include "PlcMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

class TestMappableGridGenerator : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestMakeGeometry() throw(Exception)
    {
        MappableGridGenerator<3> generator;
        boost::shared_ptr<Part<3> > p_part = generator.GenerateHemisphere(1.5, 0.1 , M_PI, M_PI/2.0, 10, 10);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part);
        VtkMeshWriter<3, 3> mesh_writer("TestMappableGridGenerator", "Hemisphere", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTMAPPABLEGRIDGENERATORHPP_*/
