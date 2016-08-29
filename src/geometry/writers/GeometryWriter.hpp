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

#ifndef GeometryWriter_HPP_
#define GeometryWriter_HPP_

#include <string>
#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK

/**
 * This class writes vtk polydata to file in VTK XML or ASCII STL format
 */
class GeometryWriter
{

private:

    /**
     * A vtk representation of the network
     */
    vtkSmartPointer<vtkPolyData> mpInputSurface;

    /**
     * The output file name
     */
    std::string mFilename;

    bool mWriteStl;

public:

    /**
     * Constructor
     */
    GeometryWriter();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     * @return a shared pointer to the class
     */
    static boost::shared_ptr<GeometryWriter> Create();

    /**
     * Destructor
     */
    ~GeometryWriter();

    /**
     * Set the polydata to be written
     */
    void SetInput(vtkSmartPointer<vtkPolyData> pSurface);

    /**
     * Adds a collection of vessels to the VesselNetwork
     * @param rFileName the full output path
     */
    void SetFileName(const std::string& rFileName);

    void SetWriteStl(bool writeStl);

    /**
     * Do the write
     */
    void Write();

};

#endif /* GeometryWriter_HPP_ */
