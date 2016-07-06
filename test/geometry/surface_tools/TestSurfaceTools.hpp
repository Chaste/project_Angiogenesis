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

#ifndef TestImageToSurface_HPP_
#define TestImageToSurface_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
//#ifdef CHASTE_ANGIOGENESIS_VMTK
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include "VtkSurfaceWriter.hpp"
#include "VtkBoundaryExtractor.hpp"
//#endif /*CHASTE_ANGIOGENESIS_VMTK*/

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"

class TestSurfaceTools : public CxxTest::TestSuite
{
public:

    void TestExtractBoundary()
    {
//        #ifdef CHASTE_ANGIOGENESIS_VMTK

        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestSurfaceTools/");
        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/surface.vtp", RelativeTo::ChasteSourceRoot);

        vtkSmartPointer<vtkXMLPolyDataReader> p_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        p_reader->SetFileName(finder.GetAbsolutePath().c_str());
        p_reader->Update();

        VtkBoundaryExtractor extractor;
        extractor.SetInput(p_reader->GetOutput());
        extractor.SetDoSmoothing(false);
        extractor.Update();

        boost::shared_ptr<VtkSurfaceWriter> p_writer = VtkSurfaceWriter::Create();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary.vtp").c_str());
        p_writer->SetInput(extractor.GetOutput());
        p_writer->Write();

        extractor.SetDoSmoothing(true);
        extractor.SetSmoothingLength(200.0);
        extractor.Update();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary_smoothed.vtp").c_str());
        p_writer->SetInput(extractor.GetOutput());
        p_writer->Write();
//        #endif /*CHASTE_ANGIOGENESIS_VMTK*/
    }
};
#endif