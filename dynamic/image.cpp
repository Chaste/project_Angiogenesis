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

#ifdef CHASTE_ANGIOGENESIS_PYTHON
#include <vector>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include "SmartPointers.hpp"
#include "ImageIO.hpp"
#include "ImageToSurface.hpp"
#include "NetworkToImage.hpp"
using namespace boost::python;

// Make the module
BOOST_PYTHON_MODULE(_image)
{
    class_<ImageIO, boost::shared_ptr<ImageIO> >("ImageIO")
        .def("SetFilename",&ImageIO::SetFilename)
        .def("SetImageResizeFactors",&ImageIO::SetImageResizeFactors)
        .def("Read",&ImageIO::Read)
        .def("Write",&ImageIO::Write)
        .def("GetImage",&ImageIO::GetImage)
        .def("SetImage",&ImageIO::SetImage)
    ;

    class_<ImageToSurface, boost::shared_ptr<ImageToSurface> >("ImageToSurface")
        .def("SetInput",&ImageToSurface::SetInputRaw)
        .def("SetThreshold",&ImageToSurface::SetThreshold)
        .def("SetUseMarchingCubes",&ImageToSurface::SetUseMarchingCubes)
        .def("Update",&ImageToSurface::Update)
        .def("GetOutput",&ImageToSurface::GetOutput)
    ;

    class_<NetworkToImage<3>, boost::shared_ptr<NetworkToImage<3> > >("NetworkToImage")
        .def("SetNetwork",&NetworkToImage<3>::SetNetwork)
        .def("SetGridSpacing",&NetworkToImage<3>::SetGridSpacing)
        .def("SetPaddingFactors",&NetworkToImage<3>::SetPaddingFactors)
        .def("SetImageDimension",&NetworkToImage<3>::SetImageDimension)
        .def("Update",&NetworkToImage<3>::Update)
        .def("GetOutput",&NetworkToImage<3>::GetOutput)
    ;
}

#endif // CHASTE_ANGIOGENESIS_PYTHON
