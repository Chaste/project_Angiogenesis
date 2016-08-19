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
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "SmartPointers.hpp"
#include "PythonObjectConverters.hpp"
#include "Vertex.hpp"
#include "Polygon.hpp"
#include "Facet.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "Cell.hpp"
#include "SimpleCell.hpp"

namespace bp = boost::python;

/**
 *  This is the preload angiogenesis module. It contains functionality that should be loaded/registered when the Chaste Python
 *  Angiogenesis package is first imported. This includes a collection of C++-Python converters.
 */

// Make the module
BOOST_PYTHON_MODULE(_chaste_project_Angiogenesis_preload)
{
    // Iterators
    PythonIterableToStl()
      .from_python<std::vector<boost::shared_ptr<Vertex> > > ()
      .from_python<std::vector<boost::shared_ptr<Polygon> > > ()
      .from_python<std::vector<boost::shared_ptr<Facet> > > ()
      .from_python<std::vector<boost::shared_ptr<VesselNode<3> > > > ()
      .from_python<std::vector<boost::shared_ptr<VesselSegment<3> > > >()
      .from_python<std::vector<boost::shared_ptr<Vessel<3> > > >()
      .from_python<std::vector<boost::shared_ptr<Cell > > >()
      .from_python<std::vector<boost::shared_ptr<SimpleCell<3> > > >()
      ;

    // Odds and ends
    bp::class_<std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > > > ("NodePair")
    .def_readwrite("first", &std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > >::first)
    .def_readwrite("second", &std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > >::second)
    ;

    bp::class_<std::pair<boost::shared_ptr<VesselSegment<3> >, double > > ("SegDistPair")
    .def_readwrite("first", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::first)
    .def_readwrite("second", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::second)
    ;
}
