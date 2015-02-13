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

#ifndef VASCULARNETWORKGENERATOR_HPP_
#define VASCULARNETWORKGENERATOR_HPP_

#include <vector>
#include <string>
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "VascularNode.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK

template<unsigned DIM>
class VasculatureGenerator
{

public:

    /*
     * Constructor
     *
     * @param prototype prototype vessel using which other vessel objects will be instantiated.
     */
    VasculatureGenerator();

    /*
     * Destructor
     */
    ~VasculatureGenerator();

    /*
     * Pattern Unit
     */
    void PatternUnitByTranslation(boost::shared_ptr<CaVascularNetwork<DIM> > input_unit,
    																		unsigned doublings_in_direction1,
    																		unsigned doublings_in_direction2,
    																		unsigned doublings_in_direction3 = 1);

    /*
     * Creates a hexagonal repeating unit
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > GenerateHexagonalUnit(double vessel_length);

    /*
     * Creates a bifurcation repeating unit
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > GenerateBifurcationUnit(double vessel_length);


	#ifdef CHASTE_VTK
    /*
	 * Generates a vessel network from a vtk file.
	 *
	 * @param filename name of file in which vascular network is described.
	 * @return a pointer to the generated vascular network.
	 */
	boost::shared_ptr<CaVascularNetwork<DIM> > GenerateNetworkFromVtkFile(const std::string& filename);
	#endif // CHASTE_VTK

};

#endif /* VASCULARNETWORKGENERATOR_HPP_ */
