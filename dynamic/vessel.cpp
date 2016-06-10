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
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VasculatureData.hpp"
#include "VascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SegmentFlowProperties.hpp"
#include "NodeFlowProperties.hpp"
#include "FlowSolver.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#include "converters.hpp"

using namespace boost::python;

// Vascular Node Overloads
boost::shared_ptr<VascularNode<3> > (*VN3_ChastePoint)(const ChastePoint<3>& ) = &VascularNode<3>::Create;
boost::shared_ptr<VascularNode<3> > (*VN3_Doubles)(double, double, double) = &VascularNode<3>::Create;
boost::shared_ptr<VascularNode<3> > (*VN3_CVec)(c_vector<double, 3>) = &VascularNode<3>::Create;
void (VascularNode<3>::*VN3_SetLocation)(c_vector<double, 3>) = &VascularNode<3>::SetLocation;

// Vessel Segment Overloads
boost::shared_ptr<VesselSegment<3> > (*VS3_Nodes)(boost::shared_ptr<VascularNode<3> >, boost::shared_ptr<VascularNode<3> >) = &VesselSegment<3>::Create;
boost::shared_ptr<VesselSegment<3> > (*VS3_Copy)(boost::shared_ptr<VesselSegment<3> >) = &VesselSegment<3>::Create;
double (VesselSegment<3>::*VS3_GetDistance)(c_vector<double, 3>) const = &VesselSegment<3>::GetDistance;
//c_vector<double, 3> (VesselSegment<3>::*VS3_GetPointProjection)(c_vector<double, 3>) const = &VesselSegment<3>::GetPointProjection;

// Vessel Overloads
boost::shared_ptr<Vessel<3> > (*V3_SingleSegment)(boost::shared_ptr<VesselSegment<3> >) = &Vessel<3>::Create;
boost::shared_ptr<Vessel<3> > (*V3_MultiSegment)(std::vector<boost::shared_ptr<VesselSegment<3> > >) = &Vessel<3>::Create;
boost::shared_ptr<Vessel<3> > (*V3_Nodes)(std::vector<boost::shared_ptr<VascularNode<3> > >) = &Vessel<3>::Create;

// Vessel Network Overloads
boost::shared_ptr<VascularNode<3> > (VascularNetwork<3>::*VNet3_NearestNode)(c_vector<double, 3>) = &VascularNetwork<3>::GetNearestNode;
std::pair<boost::shared_ptr<VesselSegment<3> >, double> (VascularNetwork<3>::*VNet3_NearestSegment)(c_vector<double, 3>) = &VascularNetwork<3>::GetNearestSegment;
boost::shared_ptr<Vessel<3> > (VascularNetwork<3>::*VNet3_NearestVessel)(c_vector<double, 3>) = &VascularNetwork<3>::GetNearestVessel;
void (VascularNetwork<3>::*VNet3_MergeCoincident)(double) = &VascularNetwork<3>::MergeCoincidentNodes;
void (VascularNetwork<3>::*VNet3_Translate)(const c_vector<double, 3>&) = &VascularNetwork<3>::Translate;

// Make the module
BOOST_PYTHON_MODULE(_vessel)
{
    class_<VascularNode<3>, boost::shared_ptr<VascularNode<3> >, boost::noncopyable >("VascularNode", no_init)
        .def("__init__", make_constructor(VN3_ChastePoint))
        .def("__init__", make_constructor(VN3_Doubles))
        .def("__init__", make_constructor(VN3_CVec))
        .add_property("id", &VascularNode<3>::GetId, &VascularNode<3>::SetId)
        .def("GetNumberOfSegments", &VascularNode<3>::GetNumberOfSegments)
        .def("GetRadius", &VascularNode<3>::GetRadius)
        .def("SetRadius", &VascularNode<3>::SetRadius)
        .def("GetVesselSegment", &VascularNode<3>::GetVesselSegment)
        .def("GetVesselSegments", &VascularNode<3>::GetVesselSegments)
        .def("GetLocation", &VascularNode<3>::GetLocationVector)
        .def("SetLocation", VN3_SetLocation)
        .def("GetFlowProperties", &VascularNode<3>::GetFlowProperties)
    ;

    class_<NodeFlowProperties, boost::shared_ptr<NodeFlowProperties>, boost::noncopyable >("NodeFlowProperties")
        .def("GetPressure", &NodeFlowProperties::GetPressure)
        .def("SetPressure", &NodeFlowProperties::SetPressure)
        .def("SetIsInputNode", &NodeFlowProperties::SetIsInputNode)
        .def("SetIsOutputNode", &NodeFlowProperties::SetIsOutputNode)
        .def("IsInputNode", &NodeFlowProperties::IsInputNode)
        .def("IsOutputNode", &NodeFlowProperties::IsOutputNode)
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

    class_<SegmentFlowProperties, boost::shared_ptr<SegmentFlowProperties>, boost::noncopyable >("SegmentFlowProperties")
        .def("SetImpedance", &SegmentFlowProperties::SetImpedance)
        .def("SetViscosity", &SegmentFlowProperties::SetViscosity)
        .def("GetFlowRate", &SegmentFlowProperties::GetFlowRate)
//        .def("GetFlowVelocity", &SegmentFlowProperties::GetFlowVelocity)
        .def("SetHaematocrit", &SegmentFlowProperties::SetHaematocrit)
        .def("GetHaematocrit", &SegmentFlowProperties::GetHaematocrit)
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
        .def("GetSegment", &Vessel<3>::GetSegment)
        .def("GetDistance", &Vessel<3>::GetDistance)
        .def("GetFlowRate", &Vessel<3>::GetFlowRate)
        .def("IsConnected", &Vessel<3>::IsConnectedTo)
        .def("RemoveSegments", &Vessel<3>::RemoveSegments)
    ;

    class_<VascularNetwork<3>, boost::shared_ptr<VascularNetwork<3> >, boost::noncopyable >("VascularNetwork")
        .def("GetNodes", &VascularNetwork<3>::GetNodes)
		.def("GetVesselEndNodes", &VascularNetwork<3>::GetVesselEndNodes)
		.def("GetVessels", &VascularNetwork<3>::GetVessels)
		.def("GetExtents", &VascularNetwork<3>::GetExtents)
		.def("GetNumberOfNodes", &VascularNetwork<3>::GetNumberOfNodes)
		.def("GetNumberOfVessels", &VascularNetwork<3>::GetNumberOfVessels)
		.def("AddVessel", &VascularNetwork<3>::AddVessel)
        .def("AddVessels", &VascularNetwork<3>::AddVessels)
        .def("CopySegmentFlowProperties", &VascularNetwork<3>::CopySegmentFlowProperties)
        .def("GetVessel", &VascularNetwork<3>::GetVessel)
        .def("GetNearestNode", VNet3_NearestNode)
        .def("GetNearestSegment", VNet3_NearestSegment)
        .def("GetNearestVessel", VNet3_NearestVessel)
        .def("MergeNodes", VNet3_MergeCoincident)
        .def("SetNodeData", &VascularNetwork<3>::SetNodeData)
        .def("GetAverageVesselLength", &VascularNetwork<3>::GetAverageVesselLength)
        .def("SetVesselData", &VascularNetwork<3>::SetVesselData)
        .def("SetSegmentRadii", &VascularNetwork<3>::SetSegmentRadii)
        .def("SetNodeRadii", &VascularNetwork<3>::SetNodeRadii)
        .def("RemoveVessel", &VascularNetwork<3>::RemoveVessel)
        .def("RemoveShortVessels", &VascularNetwork<3>::RemoveShortVessels)
        .def("MergeShortVessels", &VascularNetwork<3>::MergeShortVessels)
        .def("UpdateAll", &VascularNetwork<3>::UpdateAll)
        .def("UpdateNodes", &VascularNetwork<3>::UpdateNodes)
        .def("Translate", VNet3_Translate)
        .def("Write", &VascularNetwork<3>::Write)
        .def("GetVtk", &VascularNetwork<3>::GetVtk)
        .def("GetTotalLength", &VascularNetwork<3>::GetTotalLength)
        .def("WriteConnectivity", &VascularNetwork<3>::WriteConnectivity)
        .def("GetVesselLengthDistribution", &VascularNetwork<3>::GetVesselLengthDistribution)
        .enable_pickling()
    ;

    class_<VasculatureData>("VasculatureData")
		.def("GetScalarData", &VasculatureData::GetData<double>)
		.def("GetVectorData", &VasculatureData::GetData<std::vector<double> >)
		.def("SetScalarData", &VasculatureData::SetData<double>)
		.def("SetVectorData", &VasculatureData::SetData<std::vector<double> >)
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

    class_<std::pair<boost::shared_ptr<VascularNode<3> >, boost::shared_ptr<VascularNode<3> > > > ("NodePair")
		.def_readwrite("first", &std::pair<boost::shared_ptr<VascularNode<3> >, boost::shared_ptr<VascularNode<3> > >::first)
    	.def_readwrite("second", &std::pair<boost::shared_ptr<VascularNode<3> >, boost::shared_ptr<VascularNode<3> > >::second)
    ;

    class_<std::pair<boost::shared_ptr<VesselSegment<3> >, double > > ("SegDistPair")
		.def_readwrite("first", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::first)
    	.def_readwrite("second", &std::pair<boost::shared_ptr<VesselSegment<3> >, double >::second)
    ;

    enum_<SegmentLocation::Value>("SegmentLocation")
        .value("start", SegmentLocation::Start)
        .value("end", SegmentLocation::End)
        ;

    class_<std::vector<boost::shared_ptr<VascularNode<3> > > > ("VecNodePtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<VascularNode<3> > > >())
    ;

    class_<std::vector<boost::shared_ptr<VesselSegment<3> > > > ("VecSegmentPtrs")
        .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<VesselSegment<3> > > >())
    ;

    class_<std::vector<boost::shared_ptr<Vessel<3> > > > ("VecVesselPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Vessel<3> > > >())
    ;

    // Register Angiogenesis Converters
    PythonIterableToStl()
        .from_python<std::vector<boost::shared_ptr<VascularNode<3> > > > ()
        .from_python<std::vector<boost::shared_ptr<VesselSegment<3> > > >()
        .from_python<std::vector<boost::shared_ptr<Vessel<3> > > >()
      ;
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
