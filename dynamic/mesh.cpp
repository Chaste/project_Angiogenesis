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
#include "HybridMesh.hpp"
#include "PottsMesh.hpp"
#include "SharedPottsMeshGenerator.hpp"
#include "PottsElement.hpp"
#include "RegularGrid.hpp"
#include "Node.hpp"
#include "DistanceMap.hpp"
#include "FunctionMap.hpp"
#include "AbstractRegularGridHybridSolver.hpp"
#include "AbstractHybridSolver.hpp"

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(HyrbidMeshGenerateFromPartOverloads, HybridMesh<3>::GenerateFromPart, 1, 2);

// Make the module
BOOST_PYTHON_MODULE(_mesh)
{
    class_<RegularGrid<3>, boost::shared_ptr<RegularGrid<3> > >("RegularGrid")
        .def("GenerateFromPart", &RegularGrid<3>::GenerateFromPart)
        .def("GetExtents", &RegularGrid<3>::GetExtents)
        .def("SetExtents", &RegularGrid<3>::SetExtents)
        .def("SetSpacing", &RegularGrid<3>::SetSpacing)
        .def("SetOrigin", &RegularGrid<3>::SetOrigin)
        .def("GetSpacing", &RegularGrid<3>::GetSpacing)
        .def("GetLocations", &RegularGrid<3>::GetLocations)
        .def("GetLocationOf1dIndex", &RegularGrid<3>::GetLocationOf1dIndex)
        .def("GetNumberOfPoints", &RegularGrid<3>::GetNumberOfPoints)
        .def("GetNearestGridIndex", &RegularGrid<3>::GetNearestGridIndex)
        .def("SetPointValues", &RegularGrid<3>::SetPointValues)
        .def("SetVesselNetwork", &RegularGrid<3>::SetVesselNetwork)
        .def("SetCellPopulation", &RegularGrid<3>::SetCellPopulation)
        .def("SetCaBasedPopulation", &RegularGrid<3>::SetCaBasedPopulation)
        .def("Write", &RegularGrid<3>::Write)
        .def("SetUpVtkGrid", &RegularGrid<3>::SetUpVtkGrid)
        .def("GetVtkGrid", &RegularGrid<3>::GetVtkGrid);

    class_<HybridMesh<3>, boost::shared_ptr<HybridMesh<3> >, boost::noncopyable>("HybridMesh")
        .def("GetConnectivity", &HybridMesh<3>::GetConnectivity)
        .def("GetNodeLocations", &HybridMesh<3>::GetNodeLocations)
        .def("GenerateFromPart", &HybridMesh<3>::GenerateFromPart, HyrbidMeshGenerateFromPartOverloads());

    class_<SharedPottsMeshGenerator<3>, boost::shared_ptr<SharedPottsMeshGenerator<3> > >("PottsMeshGenerator",
                                                                              init<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned,
                                                                              unsigned, unsigned, unsigned, bool, bool, bool, bool>())
             .def("GetMesh", &SharedPottsMeshGenerator<3>::GetMesh);

    class_<PottsMesh<3u>, boost::shared_ptr<PottsMesh<3u> >, boost::noncopyable>("PottsMesh",
            init<std::vector<Node<3>*>, std::vector<PottsElement<3>*>, std::vector<std::set<unsigned> >, std::vector<std::set<unsigned> > >())
            .def("GetNumNodes", &PottsMesh<3u>::GetNumNodes)
            .def("GetNumElements", &PottsMesh<3u>::GetNumElements);
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
