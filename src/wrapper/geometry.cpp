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

#include <vector>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "SmartPointers.hpp"
#include "Vertex.hpp"
#include "Polygon.hpp"
#include "Facet.hpp"
#include "Part.hpp"
#include "MappableGridGenerator.hpp"
#include "converters.hpp"
using namespace boost::python;

// Pointer to overloaded member functions
const c_vector<double, 3>& (Vertex::*Vertex_rGetLocation)() const = &Vertex::rGetLocation;

// Autogeneration of thin wrappers for overloaded member functions
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PlcPartAddCuboidOverLoads, Part<3>::AddCuboid, 1, 4);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PlcPartAddCircleOverLoads, Part<3>::AddCircle, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PlcPartAddRectangleOverLoads, Part<3>::AddRectangle, 1, 3);
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PlcPartAddPolygonOverLoads, Part::AddPolygon, 1, 3);

//boost::shared_ptr<Polygon> f1(std::vector<boost::shared_ptr<Vertex> > x) { return &Part::AddPolygon(x); }
//boost::shared_ptr<Polygon>(Part::*AddPolygonVertices)(std::vector<boost::shared_ptr<Vertex> >, bool) = &Part::AddPolygon;

// Make the module
BOOST_PYTHON_MODULE(_geometry)
{
    class_<Vertex, boost::shared_ptr<Vertex > >("Vertex", init<optional<double, double, double> >())
        .def(init<c_vector<double, 3> >())
        .add_property("id", &Vertex::GetIndex, &Vertex::SetIndex)
        .def("GetLocation", Vertex_rGetLocation, return_value_policy<copy_const_reference>())
        .def("Translate", &Vertex::Translate)
        .def("RotateAboutAxis", &Vertex::RotateAboutAxis)
    ;

    class_<Polygon, boost::shared_ptr<Polygon> >("Polygon", init<std::vector<boost::shared_ptr<Vertex> > >())
        .def(init<boost::shared_ptr<Vertex> >())
        .def("AddVertices", &Polygon::AddVertices)
        .def("AddVertex", &Polygon::AddVertex)
        .def("GetVertices", &Polygon::GetVertices)
        .def("GetCentroid", &Polygon::GetCentroid)
        .def("GetNormal", &Polygon::GetNormal)
        .def("GetBoundingBox", &Polygon::GetBoundingBox)
        .def("GetDistanceToEdges", &Polygon::GetDistanceToEdges)
        .def("ContainsPoint", &Polygon::ContainsPoint)
        .def("Translate", &Polygon::Translate)
        .def("RotateAboutAxis", &Polygon::RotateAboutAxis)
    ;

    class_<Facet, boost::shared_ptr<Facet> >("Facet", init<std::vector<boost::shared_ptr<Polygon> > >())
        .def(init<boost::shared_ptr<Polygon> >())
        .def("AddPolygons", &Facet::AddPolygons)
        .def("AddPolygon", &Facet::AddPolygon)
        .def("GetPolygons", &Facet::GetPolygons)
        .def("GetVertices", &Facet::GetVertices)
        .def("GetCentroid", &Facet::GetCentroid)
        .def("GetNormal", &Facet::GetNormal)
        .def("UpdateVertices", &Facet::UpdateVertices)
        .def("Translate", &Facet::Translate)
        .def("RotateAboutAxis", &Facet::RotateAboutAxis)
        .def("SetData", &Facet::SetData)
    ;

    class_<Part<3>, boost::shared_ptr<Part<3> > >("Part")
        .def("AddCircle", &Part<3>::AddCircle, PlcPartAddCircleOverLoads())
        .def("AddCuboid", &Part<3>::AddCuboid, PlcPartAddCuboidOverLoads())
        //.def("add_polygon", f1)
        .def("AddRectangle", &Part<3>::AddRectangle, PlcPartAddRectangleOverLoads())
        .def("AddVesselNetwork", &Part<3>::AddVesselNetwork)
        .def("Extrude", &Part<3>::Extrude)
        .def("GetVertices", &Part<3>::GetVertices)
        .def("GetPolygons", &Part<3>::GetPolygons)
        .def("GetFacets", &Part<3>::GetFacets)
        .def("GetBoundingBox", &Part<3>::GetBoundingBox)
        .def("IsPointInPart", &Part<3>::IsPointInPart)
        .def("GetVtk", &Part<3>::GetVtk)
        .def("Translate", &Part<3>::Translate)
        .def("AddCylinder", &Part<3>::AddCylinder)
        .def("Write", &Part<3>::Write)
    ;

    class_<MappableGridGenerator>("MappableGridGenerator")
        .def("GeneratePlane", &MappableGridGenerator::GeneratePlane)
        .def("GenerateCylinder", &MappableGridGenerator::GenerateCylinder)
        .def("GenerateHemisphere", &MappableGridGenerator::GenerateHemisphere)
    ;

    // Containers
    class_<std::vector<boost::shared_ptr<Vertex> > > ("VecVertexPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Vertex> > >())
    ;
    class_<std::vector<boost::shared_ptr<Polygon> > > ("VecPolygonPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Polygon> > >())
    ;
    class_<std::vector<boost::shared_ptr<Facet> > > ("VecFacetPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Facet> > >())
    ;

    PythonIterableToStl()
      .from_python<std::vector<boost::shared_ptr<Vertex> > > ()
      .from_python<std::vector<boost::shared_ptr<Polygon> > > ()
      .from_python<std::vector<boost::shared_ptr<Facet> > > ()
      ;
}
