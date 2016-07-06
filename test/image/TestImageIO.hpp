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
#include "ImageIO.hpp"
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "Debug.hpp"

class TestImageIO : public CxxTest::TestSuite
{
public:

    void TestDefaultExtraction()
    {
        MARK;
        #ifdef CHASTE_ANGIOGENESIS_EXTENDED

        MARK;
        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestImageIO/");
        FileFinder finder = FileFinder("projects/Angiogenesis/test/data/median.tif", RelativeTo::ChasteSourceRoot);

        std::cout << "wut!!" << std::endl;
        MARK;
        ImageIO image_io;

        MARK;
        // Read the file in vtk format
        image_io.SetFilename(finder.GetAbsolutePath());
        MARK;
        image_io.SetImageResizeFactors(0.5, 0.5, 1.0);
        MARK;
        image_io.ReadVtkImage();

        MARK;
        // Write it out in VTI format
        image_io.SetFilename(file_handler1.GetOutputDirectoryFullPath()+"image_vtk_format.vti");
        image_io.WriteVtkImage();

        MARK;
        // Write it out in PNG format using ITK
        image_io.SetFilename(file_handler1.GetOutputDirectoryFullPath()+"image_itk_format.png");
        MARK;
        image_io.ConvertVtkToItk();
        MARK;
        image_io.WriteItkImage();
        MARK;

        // Read the image using ITK
        image_io.SetFilename(finder.GetAbsolutePath());
        MARK;
        image_io.ReadItkImage();
        MARK;

        // Write it out in VTK format
        image_io.ConvertItkToVtk();
        image_io.SetFilename(file_handler1.GetOutputDirectoryFullPath()+"image_vtk_from_itk_format.vti");
        image_io.WriteVtkImage();

        #endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
    }
};
#endif