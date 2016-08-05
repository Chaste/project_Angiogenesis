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

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkVersion.h>
#endif // CHASTE_VTK
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "VtkSurfaceWriter.hpp"

VtkSurfaceWriter::VtkSurfaceWriter() :
    mpInputSurface(),
    mFilename(),
    mWriteStl(false)
{

}

VtkSurfaceWriter::~VtkSurfaceWriter()
{

}

boost::shared_ptr<VtkSurfaceWriter > VtkSurfaceWriter::Create()
{
    MAKE_PTR(VtkSurfaceWriter, pSelf);
    return pSelf;
}

void VtkSurfaceWriter::SetInput(vtkSmartPointer<vtkPolyData> pSurface)
{
    mpInputSurface = pSurface;
}

void VtkSurfaceWriter::SetWriteStl(bool writeStl)
{
    mWriteStl = writeStl;
}

void VtkSurfaceWriter::SetFileName(const std::string& rFileName)
{
    mFilename = rFileName;
}

void VtkSurfaceWriter::Write()
{
    if(!mpInputSurface)
    {
        EXCEPTION("An input surface is not set.");
    }

    if(mFilename.empty())
    {
        EXCEPTION("No file name set for VtkVesselNetworkWriter");
    }

    if(mWriteStl)
    {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
        writer->SetFileName(mFilename.c_str());
        #if VTK_MAJOR_VERSION <= 5
            writer->SetInput(mpInputSurface);
        #else
            writer->SetInputData(mpInputSurface);
        #endif
        writer->SetFileTypeToASCII();
        writer->Write();
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(mFilename.c_str());
        #if VTK_MAJOR_VERSION <= 5
            writer->SetInput(mpInputSurface);
        #else
            writer->SetInputData(mpInputSurface);
        #endif
        writer->Write();
    }
}
