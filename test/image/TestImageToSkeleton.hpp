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
#ifdef CHASTE_ANGIOGENESIS_EXTENDED
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include "ImageToSkeleton.hpp"
#include "ImageIO.hpp"
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"

class TestImageToSkeleton : public CxxTest::TestSuite
{
public:

    void TestDefaultExtraction()
    {
        #ifdef CHASTE_ANGIOGENESIS_EXTENDED

        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestImageToSkeleton/");
//        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/median.tif", RelativeTo::ChasteSourceRoot);
//        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/domain.png", RelativeTo::ChasteSourceRoot);
//        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/test.tif", RelativeTo::ChasteSourceRoot);

        ImageIO reader = ImageIO();
//        reader.SetFilename(finder.GetAbsolutePath());
        reader.SetFilename("/home/grogan/median.tif");
        reader.SetImageResizeFactors(0.5, 0.5, 1.0);
        reader.ReadVtkImage();

        vtkSmartPointer<vtkXMLImageDataWriter> p_writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"initial.vti").c_str());
        p_writer->SetInput(reader.GetVtkImage());
        p_writer->Write();

        // Extract the surface
        ImageToSkeleton skeleton_extract = ImageToSkeleton();
        skeleton_extract.SetInput(reader.GetVtkImage());
        skeleton_extract.SetReverseIntensity(true);
        skeleton_extract.SetUseVtkVersion(false);
        skeleton_extract.Update();

        vtkSmartPointer<vtkXMLImageDataWriter> p_writer1 = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        p_writer1->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"skeleton.vti").c_str());
        p_writer1->SetInput(skeleton_extract.GetOutput());
        p_writer1->Write();
        #endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
    }
};
#endif