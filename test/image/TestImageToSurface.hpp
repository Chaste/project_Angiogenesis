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
#include "ImageToSurface.hpp"
#include "ImageReader.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>

class TestImageToSurface : public CxxTest::TestSuite
{
public:

    void TestDefaultExtraction()
    {
        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestImageToSurface/");
        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/median.tif", RelativeTo::ChasteSourceRoot);

        ImageReader reader = ImageReader();
        reader.SetFilename(finder.GetAbsolutePath());
        reader.SetImageResizeFactors(0.5, 0.5, 1.0);
        reader.Update();

        vtkSmartPointer<vtkXMLImageDataWriter> p_writer1 = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        p_writer1->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"image.vti").c_str());
        p_writer1->SetInput(reader.GetOutput());
        p_writer1->Write();

        // Extract the surface
        ImageToSurface surface_extract = ImageToSurface();
        surface_extract.SetInput(reader.GetOutput());
        surface_extract.SetThreshold(1.0, false);
        surface_extract.Update();

        // Write the surface to file
        vtkSmartPointer<vtkXMLPolyDataWriter> p_writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"surface.vtp").c_str());
        p_writer->SetInput(surface_extract.GetOutput());
        p_writer->Write();
    }
};

#endif /*TestImageToSurface_HPP_*/