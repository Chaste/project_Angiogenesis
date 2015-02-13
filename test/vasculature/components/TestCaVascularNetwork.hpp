/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TESTVESSELNETWORK_HPP_
#define TESTVESSELNETWORK_HPP_

#include <math.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "VasculatureData.hpp"
#include "ChastePoint.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"
#include "OutputFileHandler.hpp"
#include "FakePetscSetup.hpp"

class TestVesselNetwork : public AbstractCellBasedTestSuite
{
public:

	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

    void TestConstructor() throw(Exception)
    {
    	// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(1.0, 2.0, 6.0));
		points.push_back(ChastePoint<3>(3.0, 4.0, 7.0));
		points.push_back(ChastePoint<3>(3.0, 4.0, 7.0));
		points.push_back(ChastePoint<3>(3.0, 4.0, 8.0));
		points.push_back(ChastePoint<3>(3.0, 4.0, 9.0));

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}
		nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[1])));

    	// Make some segments
		SegmentPtr3 pSegment1(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
		SegmentPtr3 pSegment2(CaVesselSegment<3>::Create(nodes[2], nodes[3]));
		SegmentPtr3 pSegment3(CaVesselSegment<3>::Create(nodes[3], nodes[4]));

		// Make some vessels
		VesselPtr3 pVessel1(CaVessel<3>::Create(pSegment1));
		VesselPtr3 pVessel2(CaVessel<3>::Create(pSegment2));
		VesselPtr3 pVessel3(CaVessel<3>::Create(pSegment3));

    	std::vector<VesselPtr3> vessels;
    	vessels.push_back(pVessel2);
    	vessels.push_back(pVessel3);

		// Make a network
		CaVascularNetwork<3> vessel_network;
		vessel_network.AddVessel(pVessel1);
		vessel_network.AddVessels(vessels);

		TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 5u);

		vessel_network.MergeCoincidentNodes();

		TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 4u);

		// Try writing to file
		OutputFileHandler output_file_handler("TestVesselNetwork", false);
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("GenericVesselNetwork.vtp");
		vessel_network.WriteToFile(output_filename, true);

		// Move the network
		std::vector<double> translation_vector_3d;
		translation_vector_3d.push_back(3.5);
		translation_vector_3d.push_back(5.6);
		translation_vector_3d.push_back(-12.8);

		vessel_network.Translate(translation_vector_3d);
		vessel_network.Translate(translation_vector_3d, true);
		std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("GenericVesselNetwork_TranslatedCopy.vtp");
		vessel_network.WriteToFile(output_filename2, true);
    }
};

#endif /*TESTVESSELNETWORK_HPP_*/
