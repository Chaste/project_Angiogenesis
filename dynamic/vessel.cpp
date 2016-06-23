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
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "UblasIncludes.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SegmentFlowProperties.hpp"
#include "NodeFlowProperties.hpp"
#include "VesselFlowProperties.hpp"
#include "converters.hpp"

using namespace boost::python;

// Vascular Node Overloads
boost::shared_ptr<VesselNode<3> > (*VN3_Doubles)(double, double, double) = &VesselNode<3>::Create;
boost::shared_ptr<VesselNode<3> > (*VN3_CVec)(const c_vector<double, 3>&) = &VesselNode<3>::Create;
void (VesselNode<3>::*VN3_SetLocation)(const c_vector<double, 3>&) = &VesselNode<3>::SetLocation;

// Vessel Segment Overloads
boost::shared_ptr<VesselSegment<3> > (*VS3_Nodes)(boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> >) = &VesselSegment<3>::Create;
boost::shared_ptr<VesselSegment<3> > (*VS3_Copy)(boost::shared_ptr<VesselSegment<3> >) = &VesselSegment<3>::Create;
double (VesselSegment<3>::*VS3_GetDistance)(const c_vector<double, 3>&) const = &VesselSegment<3>::GetDistance;
//c_vector<double, 3> (VesselSegment<3>::*VS3_GetPointProjection)(c_vector<double, 3>) const = &VesselSegment<3>::GetPointProjection;

// Vessel Overloads
boost::shared_ptr<Vessel<3> > (*V3_SingleSegment)(boost::shared_ptr<VesselSegment<3> >) = &Vessel<3>::Create;
boost::shared_ptr<Vessel<3> > (*V3_MultiSegment)(std::vector<boost::shared_ptr<VesselSegment<3> > >) = &Vessel<3>::Create;
boost::shared_ptr<Vessel<3> > (*V3_Nodes)(std::vector<boost::shared_ptr<VesselNode<3> > >) = &Vessel<3>::Create;

// Vessel Network Overloads
boost::shared_ptr<VesselNode<3> > (VesselNetwork<3>::*VNet3_NearestNode)(const c_vector<double, 3>&) = &VesselNetwork<3>::GetNearestNode;
std::pair<boost::shared_ptr<VesselSegment<3> >, double> (VesselNetwork<3>::*VNet3_NearestSegment)(const c_vector<double, 3>&) = &VesselNetwork<3>::GetNearestSegment;
boost::shared_ptr<Vessel<3> > (VesselNetwork<3>::*VNet3_NearestVessel)(const c_vector<double, 3>&) = &VesselNetwork<3>::GetNearestVessel;
void (VesselNetwork<3>::*VNet3_SetSegmentRadii)(double) = &VesselNetwork<3>::SetSegmentRadii;
void (VesselNetwork<3>::*VNet3_MergeCoincident)(double) = &VesselNetwork<3>::MergeCoincidentNodes;
void (VesselNetwork<3>::*VNet3_Translate)(const c_vector<double, 3>&) = &VesselNetwork<3>::Translate;

// Make the module
BOOST_PYTHON_MODULE(_vessel)
{
    class_<VesselNode<3>, boost::shared_ptr<VesselNode<3> >, boost::noncopyable >("VesselNode", no_init)
        .def("__init__", make_constructor(VN3_Doubles))
        .def("__init__", make_constructor(VN3_CVec))
        .add_property("id", &VesselNode<3>::GetId, &VesselNode<3>::SetId)
        .def("GetNumberOfSegments", &VesselNode<3>::GetNumberOfSegments)
        .def("GetRadius", &VesselNode<3>::GetRadius)
        .def("SetRadius", &VesselNode<3>::SetRadius)
//        .def("GetSegment", &VesselNode<3>::GetSegment)
        .def("GetSegments", &VesselNode<3>::GetSegments)
        .def("GetLocation", &VesselNode<3>::rGetLocation, return_value_policy<copy_const_reference>())
        .def("SetLocation", VN3_SetLocation)
        .def("GetFlowProperties", &VesselNode<3>::GetFlowProperties)
    ;

    class_<NodeFlowProperties<3>, boost::shared_ptr<NodeFlowProperties<3> >, boost::noncopyable >("NodeFlowProperties")
        .def("GetPressure", &NodeFlowProperties<3>::GetPressure)
        .def("SetPressure", &NodeFlowProperties<3>::SetPressure)
        .def("SetIsInputNode", &NodeFlowProperties<3>::SetIsInputNode)
        .def("SetIsOutputNode", &NodeFlowProperties<3>::SetIsOutputNode)
        .def("IsInputNode", &NodeFlowProperties<3>::IsInputNode)
        .def("IsOutputNode", &NodeFlowProperties<3>::IsOutputNode)
    ;

    class_<VesselSegment<3>, boost::shared_ptr<VesselSegment<3> >, boost::noncopyable >("VesselSegment", no_init)
        .def("__init__", make_constructor(VS3_Copy))
        .def("__init__", make_constructor(VS3_Nodes))
        .add_property("id", &VesselSegment<3>::GetId, &VesselSegment<3>::SetId)
        .def("GetLength", &VesselSegment<3>::GetLength)
        .def("GetMidPoint", &VesselSegment<3>::GetMidPoint)
        .def("GetNodes", &VesselSegment<3>::GetNodes)
        .def("GetRadius", &VesselSegment<3>::GetRadius)
        .def("SetRadius", &VesselSegment<3>::SetRadius)
        .def("GetUnitTangent", &VesselSegment<3>::GetUnitTangent)
        .def("GetVessel", &VesselSegment<3>::GetVessel)
        .def("GetDistance", VS3_GetDistance)
        .def("GetFlowProperties", &VesselSegment<3>::GetFlowProperties)
        .def("GetOppositeNode", &VesselSegment<3>::GetOppositeNode)
        .def("ReplaceNode", &VesselSegment<3>::ReplaceNode)
    ;

    class_<SegmentFlowProperties<3>, boost::shared_ptr<SegmentFlowProperties<3> >, boost::noncopyable >("SegmentFlowProperties")
        .def("SetImpedance", &SegmentFlowProperties<3>::SetImpedance)
        .def("SetViscosity", &SegmentFlowProperties<3>::SetViscosity)
        .def("GetFlowRate", &SegmentFlowProperties<3>::GetFlowRate)
        .def("SetHaematocrit", &SegmentFlowProperties<3>::SetHaematocrit)
        .def("GetHaematocrit", &SegmentFlowProperties<3>::GetHaematocrit)
    ;

    class_<Vessel<3>, boost::shared_ptr<Vessel<3> >, boost::noncopyable >("Vessel", no_init)
        .def("__init__", make_constructor(V3_SingleSegment))
        .def("__init__", make_constructor(V3_MultiSegment))
        .def("__init__", make_constructor(V3_Nodes))
        .def("GetEndNode", &Vessel<3>::GetEndNode)
        .add_property("id", &Vessel<3>::GetId, &Vessel<3>::SetId)
        .def("GetLength", &Vessel<3>::GetLength)
        .def("GetNodes", &Vessel<3>::GetNodes)
        .def("GetNumberOfNodes", &Vessel<3>::GetNumberOfNodes)
        .def("GetNumberOfSegments", &Vessel<3>::GetNumberOfSegments)
        .def("GetSegments", &Vessel<3>::GetSegments)
        .def("GetStartNode", &Vessel<3>::GetStartNode)
        .def("AddSegment", &Vessel<3>::AddSegment)
        .def("AddSegments", &Vessel<3>::AddSegments)
//        .def("GetSegment", &Vessel<3>::GetSegment)
        .def("GetDistance", &Vessel<3>::GetDistance)
        .def("IsConnected", &Vessel<3>::IsConnectedTo)
        .def("RemoveSegments", &Vessel<3>::RemoveSegments)
    ;

    class_<VesselFlowProperties<3>, boost::shared_ptr<VesselFlowProperties<3> >, boost::noncopyable >("VesselFlowProperties")
        .def("SetImpedance", &VesselFlowProperties<3>::SetImpedance)
        .def("SetViscosity", &VesselFlowProperties<3>::SetViscosity)
        .def("GetFlowRate", &VesselFlowProperties<3>::GetFlowRate)
        .def("SetHaematocrit", &VesselFlowProperties<3>::SetHaematocrit)
        .def("GetHaematocrit", &VesselFlowProperties<3>::GetHaematocrit)
    ;

    class_<VesselNetwork<3>, boost::shared_ptr<VesselNetwork<3> >, boost::noncopyable >("VesselNetwork")
        .def("GetNodes", &VesselNetwork<3>::GetNodes)
		.def("GetVesselEndNodes", &VesselNetwork<3>::GetVesselEndNodes)
		.def("GetVessels", &VesselNetwork<3>::GetVessels)
		.def("GetExtents", &VesselNetwork<3>::GetExtents)
		.def("GetNumberOfNodes", &VesselNetwork<3>::GetNumberOfNodes)
		.def("GetNumberOfVessels", &VesselNetwork<3>::GetNumberOfVessels)
		.def("AddVessel", &VesselNetwork<3>::AddVessel)
        .def("AddVessels", &VesselNetwork<3>::AddVessels)
        .def("CopySegmentFlowProperties", &VesselNetwork<3>::CopySegmentFlowProperties)
//        .def("GetVessel", &VesselNetwork<3>::GetVessel)
        .def("GetNearestNode", VNet3_NearestNode)
        .def("GetNearestSegment", VNet3_NearestSegment)
        .def("GetNearestVessel", VNet3_NearestVessel)
        .def("MergeNodes", VNet3_MergeCoincident)
        .def("GetAverageVesselLength", &VesselNetwork<3>::GetAverageVesselLength)
        .def("SetNodeRadii", &VesselNetwork<3>::SetNodeRadii)
        .def("RemoveVessel", &VesselNetwork<3>::RemoveVessel)
        .def("RemoveShortVessels", &VesselNetwork<3>::RemoveShortVessels)
        .def("MergeShortVessels", &VesselNetwork<3>::MergeShortVessels)
        .def("UpdateAll", &VesselNetwork<3>::UpdateAll)
        .def("UpdateNodes", &VesselNetwork<3>::UpdateNodes)
        .def("SetSegmentRadii", VNet3_SetSegmentRadii)
        .def("Translate", VNet3_Translate)
        .def("Write", &VesselNetwork<3>::Write)
        .def("GetVtk", &VesselNetwork<3>::GetVtk)
        .def("GetTotalLength", &VesselNetwork<3>::GetTotalLength)
        .def("WriteConnectivity", &VesselNetwork<3>::WriteConnectivity)
        .def("GetVesselLengthDistribution", &VesselNetwork<3>::GetVesselLengthDistribution)
        .enable_pickling()
    ;

    class_<VasculatureGenerator<3> >("VasculatureGenerator")
        .def("GenerateHexagonalUnit", &VasculatureGenerator<3>::GenerateHexagonalUnit)
        .def("GenerateHexagonalNetwork", &VasculatureGenerator<3>::GenerateHexagonalNetwork)
        .def("GenerateSingleVessel", &VasculatureGenerator<3>::GenerateSingleVessel)
        .def("PatternUnitByTranslation", &VasculatureGenerator<3>::PatternUnitByTranslation)
        .def("generateBifurcationUnit", &VasculatureGenerator<3>::GenerateBifurcationUnit)
        .def("GenerateNetworkFromVtkFile", &VasculatureGenerator<3>::GenerateNetworkFromVtkFile)
        .def("GenerateParrallelNetwork", &VasculatureGenerator<3>::GenerateParrallelNetwork)
        .def("GenerateOvalNetwork", &VasculatureGenerator<3>::GenerateOvalNetwork)
    ;

    enum_<VesselDistribution::Value>("VesselDistribution")
        .value("REGULAR", VesselDistribution::REGULAR)
        .value("UNIFORM", VesselDistribution::UNIFORM)
        .value("TWO_LAYER", VesselDistribution::TWO_LAYER)
        ;

    class_<std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > > > ("NodePair")
		.def_readwrite("first", &std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > >::first)
    	.def_readwrite("second", &std::pair<boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> > >::second)
    ;

    class_<std::pair<boost::shared_ptr<VesselSegment<3> >, double > > ("SegDistPair")
		.def_readwrite("first", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::first)
    	.def_readwrite("second", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::second)
    ;

    enum_<SegmentLocation::Value>("SegmentLocation")
        .value("start", SegmentLocation::Start)
        .value("end", SegmentLocation::End)
        ;

    class_<std::vector<boost::shared_ptr<VesselNode<3> > > > ("VecNodePtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<VesselNode<3> > > >())
    ;

    class_<std::vector<boost::shared_ptr<VesselSegment<3> > > > ("VecSegmentPtrs")
        .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<VesselSegment<3> > > >())
    ;

    class_<std::vector<boost::shared_ptr<Vessel<3> > > > ("VecVesselPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Vessel<3> > > >())
    ;

    // Register Angiogenesis Converters
    PythonIterableToStl()
        .from_python<std::vector<boost::shared_ptr<VesselNode<3> > > > ()
        .from_python<std::vector<boost::shared_ptr<VesselSegment<3> > > >()
        .from_python<std::vector<boost::shared_ptr<Vessel<3> > > >()
      ;
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
