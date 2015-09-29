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

#ifndef PLCMESH_HPP_
#define PLCMESH_HPP_

#include <vector>
#include "SmartPointers.hpp"
#include "ChastePoint.hpp"
#include "TetrahedralMesh.hpp"
#include "Part.hpp"

// Jonathan Shewchuk's triangle and Hang Si's tetgen. Tetgen 1.5 is used in place of 1.4.2 used in Chaste/mesh.
// This allows more robust meshing of parts with floating segments. Note that the Tetgen 1.5 license is
// incompatible with Chaste, #2486.

#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen15.h"
#undef REAL
#undef VOID

struct triangulateio;

/**
 * A concrete TetrahedralMesh which can be constructed through
 * the input of PLC (piecewise linear complex) info.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class PlcMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
public:

    /*
     * Constructor
     */
    PlcMesh();

    /*
     * Destructor
     */
    ~PlcMesh();

    /* Factory constructor method
     * @return a shared pointer to a new mesh
     */
    static boost::shared_ptr<PlcMesh<ELEMENT_DIM, SPACE_DIM> > Create();


    void GenerateFromPart(boost::shared_ptr<Part<SPACE_DIM> > pPart, double maxElementArea = 0.0, bool useTetgen1_5 = true);

    std::vector<std::vector<unsigned> > GetConnectivity();

    std::vector<std::vector<double> > GetNodeLocations();

    void Write(const std::string& fileName);

private:

    void Mesh2d(boost::shared_ptr<Part<SPACE_DIM> > pPart, double maxElementArea = 0.0);

    void Mesh3d(boost::shared_ptr<Part<SPACE_DIM> > pPart, double maxElementArea = 0.0);

    // This is the same as the TetrahedralMesh implementation of ImportFromMesher but avoids a lot of templating hassle.
    void ImportFromTetgen(tetgen15::tetgenio& mesherOutput, unsigned numberOfElements, int *elementList,
                          unsigned numberOfFaces, int *faceList, int *edgeMarkerList);

    void InitialiseTriangulateIo(triangulateio& mesherIo);

    void FreeTriangulateIo(triangulateio& mesherIo);
};

#endif /* PLCMESH_HPP_*/
