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

#include "VasculatureGenerator.hpp"

template<unsigned DIM>
VasculatureGenerator<DIM>::VasculatureGenerator()
{
}

template<unsigned DIM>
VasculatureGenerator<DIM>::~VasculatureGenerator()
{
}

template<unsigned DIM>
void VasculatureGenerator<DIM>::PatternUnitByTranslation(boost::shared_ptr<CaVascularNetwork<DIM> > input_unit,
																		unsigned number_in_direction1,
																		unsigned number_in_direction2,
																		unsigned number_in_direction3)
{

	// Get unit dimensions
	std::set<boost::shared_ptr<VascularNode<DIM> > > nodes = input_unit->GetNodes();

	std::vector<double> max_vals;
	std::vector<double> min_vals;
	for (unsigned i = 0; i < DIM; i++)
	{
		max_vals.push_back(-1.e6);
		min_vals.push_back(1.e6);
	}

	typename std::set<boost::shared_ptr<VascularNode<DIM> > >::iterator node_iter;
	for (node_iter=nodes.begin(); node_iter!=nodes.end(); node_iter++)
	{
		for(unsigned i=0; i <DIM; i++)
		{
			if((*node_iter)->GetLocation()[i] > max_vals[i])
			{
				max_vals[i] = (*node_iter)->GetLocation()[i];
			}
			if((*node_iter)->GetLocation()[i] < min_vals[i])
			{
				min_vals[i] = (*node_iter)->GetLocation()[i];
			}
		}
	}
//
	// Get the translation vector
	std::vector<double> base_translation_vector;
	for(unsigned i=0; i <DIM; i++)
	{
		base_translation_vector.push_back(max_vals[i] - min_vals[i]);
	}

	// Copy the vessels and translate the new ones a specified number of
	// times in each direction. The translate method also merges the nodes.
	std::vector<double> translation_vector;
	translation_vector.push_back(0.0);
	translation_vector.push_back(0.0);
	if(DIM>2)
	{
		translation_vector.push_back(0.0);
	}

	for(unsigned i=0; i < number_in_direction1; i++)
	{
		if(number_in_direction1 != 1)
		{
			double distance = double(pow(2,i)) * base_translation_vector[0];
			translation_vector[0] = distance;
			input_unit->Translate(translation_vector, true);
		}
	}
	for(unsigned j=0; j < number_in_direction2; j++)
	{
		if(number_in_direction2 != 1)
		{
			double distance = double(pow(2,j)) * base_translation_vector[1];
			translation_vector[0] = 0.0;
			translation_vector[1] = distance;
			input_unit->Translate(translation_vector, true);
		}
	}
	if(DIM>2)
	{
		for(unsigned k=0; k < number_in_direction3; k++)
		{
			if(number_in_direction3 != 1)
			{
				double distance = double(pow(2,k)) * base_translation_vector[2];
				translation_vector[1] = 0.0;
				translation_vector[2] = distance;
				input_unit->Translate(translation_vector, true);
			}
		}
	}
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateHexagonalUnit(double vessel_length)
{
	// Generate the nodes
	std::vector<ChastePoint<DIM> > points;
	if(DIM == 2)
	{
		points.push_back(ChastePoint<DIM>(0.0, 0.0)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, vessel_length)); //1
		points.push_back(ChastePoint<DIM>(0.0, 2.0 * vessel_length)); //2
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, vessel_length)); //3
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length, 0.0)); //4
		points.push_back(ChastePoint<DIM>(4.0 * vessel_length, 0.0)); //5
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length,  2.0 * vessel_length)); //6
	}
	else
	{
		points.push_back(ChastePoint<DIM>(0.0, 0.0, 0.0)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, vessel_length, 0.0)); //1
		points.push_back(ChastePoint<DIM>(0.0, 2.0 * vessel_length, 0.0)); //2
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, vessel_length, 0.0)); //3
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length, 0.0, 0.0)); //4
		points.push_back(ChastePoint<DIM>(4.0 * vessel_length, 0.0, 0.0)); //5
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length,  2.0 * vessel_length, 0.0)); //6
	}

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
	typename std::vector<ChastePoint<DIM> >::iterator it;
	for(it = points.begin(); it < points.end(); it++)
	{
		nodes.push_back(boost::shared_ptr<VascularNode<DIM> >(VascularNode<DIM>::Create(*it)));
	}

	// Generate the segments and vessels
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[1], nodes[2])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[1], nodes[3])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[3], nodes[4])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[4], nodes[5])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[3], nodes[6])));

	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it2;
	for(it2 = segments.begin(); it2 < segments.end(); it2++)
	{
		vessels.push_back(boost::shared_ptr<CaVessel<DIM> > (CaVessel<DIM>::Create(*it2)));
	}

	// Generate the network
	boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork(new CaVascularNetwork<DIM>());
	pNetwork->AddVessels(vessels);

	return pNetwork;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateBifurcationUnit(double vessel_length)
{
	// Generate the nodes
	std::vector<ChastePoint<DIM> > points;
	if(DIM == 2)
	{
		points.push_back(ChastePoint<DIM>(0.0, vessel_length)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, vessel_length)); //1
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, 2.0 * vessel_length)); //2
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, 0.0)); //3
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length, vessel_length)); //4
		points.push_back(ChastePoint<DIM>(4.0 * vessel_length, vessel_length)); //5
	}
	else
	{
		points.push_back(ChastePoint<DIM>(0.0, vessel_length, 0.0)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, vessel_length, 0.0)); //1
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, 2.0 * vessel_length, 0.0)); //2
		points.push_back(ChastePoint<DIM>(2.0 * vessel_length, 0.0, 0.0)); //3
		points.push_back(ChastePoint<DIM>(3.0 * vessel_length, vessel_length, 0.0)); //4
		points.push_back(ChastePoint<DIM>(4.0 * vessel_length, vessel_length, 0.0)); //5
	}

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
	typename std::vector<ChastePoint<DIM> >::iterator it;
	for(it = points.begin(); it < points.end(); it++)
	{
		nodes.push_back((VascularNode<DIM>::Create(*it)));
	}

	// Generate the segments and vessels
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[1], nodes[2])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[1], nodes[3])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[2], nodes[4])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[3], nodes[4])));
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[4], nodes[5])));

	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it2;
	for(it2 = segments.begin(); it2 < segments.end(); it2++)
	{
		vessels.push_back(boost::shared_ptr<CaVessel<DIM> > (CaVessel<DIM>::Create(*it2)));
	}

	// Generate the network
	boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork(new CaVascularNetwork<DIM>());
	pNetwork->AddVessels(vessels);
	return pNetwork;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateSingleVessel(double vessel_length)
{
	std::vector<ChastePoint<DIM> > points;

	if(DIM == 2)
	{
		points.push_back(ChastePoint<DIM>(0.0, 0.0)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, 0.0)); //1
	}
	else
	{
		points.push_back(ChastePoint<DIM>(0.0, 0.0, 0.0)); //0
		points.push_back(ChastePoint<DIM>(vessel_length, 0.0, 0.0)); //1
	}

	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
	typename std::vector<ChastePoint<DIM> >::iterator it;
	for(it = points.begin(); it < points.end(); it++)
	{
		nodes.push_back((VascularNode<DIM>::Create(*it)));
	}

	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
	segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> >(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));

	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
	typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator it2;
	for(it2 = segments.begin(); it2 < segments.end(); it2++)
	{
		vessels.push_back(boost::shared_ptr<CaVessel<DIM> > (CaVessel<DIM>::Create(*it2)));
	}

	// Generate the network
	boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork(new CaVascularNetwork<DIM>());
	pNetwork->AddVessels(vessels);
	return pNetwork;
}

// todo This method needs to be re-implemented using new vasculature structure and naming conventions.
template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateHexagonalNetwork(unsigned width,
                unsigned height,
                unsigned vessel_length)
{

        // Determine the number of repeating hexagonal units in the x and y directions based on the input height, width and vessel length
        double horizontal_vessel_length = vessel_length;
        double diagonal_vessel_length = floor(vessel_length/std::sqrt(2));
        unsigned units_in_x_direction = floor(((width-horizontal_vessel_length-2.0*diagonal_vessel_length-1.0)/
                        (2.0*(horizontal_vessel_length+diagonal_vessel_length))));
        unsigned units_in_y_direction = floor(((height-1.0)/(2.0*diagonal_vessel_length)));

        // Ensure there are a minimal number of units required to generate a functional network
        if(units_in_x_direction <= 1 || units_in_y_direction <= 1)
        {
                EXCEPTION("Insufficient number of repeating units specified for the hexagonal network.");
        }

        // Generate an array of vessels with no spatial info
        unsigned number_of_vessels = (units_in_x_direction*units_in_y_direction*6)+
                        (units_in_y_direction*5)+units_in_x_direction;
        boost::shared_ptr<CaVascularNetwork<DIM> > pVesselNetwork(new CaVascularNetwork<DIM>());
        std::vector<boost::shared_ptr<CaVessel<DIM> > > vessel_array;
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > vessel_segment_array;

        for (unsigned i = 0; i < number_of_vessels; i++)
        {
            // instantiate set of vessel segments with any nodes then later change location of that node
            // with GetNode(index)->SetLocation(ChastePoint ... )
            boost::shared_ptr<VascularNode<DIM> > node1(VascularNode<DIM>::Create(0,0,0));
            boost::shared_ptr<VascularNode<DIM> > node2(VascularNode<DIM>::Create(0,0,0));
            vessel_segment_array.push_back(CaVesselSegment<DIM>::Create(node1, node2));
        }

        // Set the left side coordinates of vessels in the network
        vessel_segment_array[0]->GetNode(0)->SetLocation(ChastePoint<DIM>(0.0, 0.0));
        vessel_segment_array[1]->GetNode(0)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length,diagonal_vessel_length));
        vessel_segment_array[2]->GetNode(0)->SetLocation(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length, 0.0));
        for (unsigned i = 3; i < (2+3*units_in_x_direction); i++)
        {
                vessel_segment_array[i]->GetNode(0)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i-3]->GetNode(0)->GetLocation()[0] +
                                2.0*(diagonal_vessel_length+horizontal_vessel_length), vessel_segment_array[i-3]->GetNode(0)->GetLocation()[1]));
        }
        vessel_segment_array[2+3*units_in_x_direction]->GetNode(0)->SetLocation(ChastePoint<DIM>(0, 2.0*diagonal_vessel_length));
        vessel_segment_array[2+3*units_in_x_direction+1]->GetNode(0)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
        vessel_segment_array[2+3*units_in_x_direction+2]->GetNode(0)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length, diagonal_vessel_length));

        for (unsigned i = (2+3*units_in_x_direction+3); i < (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i++)
        {
                vessel_segment_array[i]->GetNode(0)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i-3]->GetNode(0)->GetLocation()[0] +
                                2.0*(diagonal_vessel_length+horizontal_vessel_length), vessel_segment_array[i-3]->GetNode(0)->GetLocation()[1]));
        }

        for (unsigned i = (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i < number_of_vessels - units_in_x_direction; i++)
        {
                vessel_segment_array[i]->GetNode(0)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode(0)->GetLocation()[0],
                                vessel_segment_array[i-(2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode(0)->GetLocation()[1] + 2.0*diagonal_vessel_length));
        }
        vessel_segment_array[number_of_vessels-units_in_x_direction]->GetNode(0)->SetLocation(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length,
                        vessel_segment_array[number_of_vessels-units_in_x_direction-1]->GetNode(0)->GetLocation()[1]+diagonal_vessel_length));

        for (unsigned i = 1; i < units_in_x_direction; i++)
        {
                vessel_segment_array[number_of_vessels-units_in_x_direction+i]->GetNode(0)->SetLocation(ChastePoint<DIM>(
                                vessel_segment_array[number_of_vessels-units_in_x_direction+i-1]->GetNode(0)->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
                                vessel_segment_array[number_of_vessels-units_in_x_direction+i-1]->GetNode(0)->GetLocation()[1]));
        }

        // Set the right side coordinates of vessels in the network
        vessel_segment_array[0]->GetNode(1)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
        vessel_segment_array[1]->GetNode(1)->SetLocation(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length, 0.0));
        vessel_segment_array[2]->GetNode(1)->SetLocation(ChastePoint<DIM>(2.0*(diagonal_vessel_length+horizontal_vessel_length), 0.0));

        for (unsigned i = 3; i < (2+3*units_in_x_direction); i++)
        {
                vessel_segment_array[i]->GetNode(1)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i-3]->GetNode(1)->GetLocation()[0]+2*(diagonal_vessel_length+horizontal_vessel_length),
                                vessel_segment_array[i-3]->GetNode(1)->GetLocation()[1]));
        }
        vessel_segment_array[2+3*units_in_x_direction]->GetNode(1)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
        vessel_segment_array[2+3*units_in_x_direction+1]->GetNode(1)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length, diagonal_vessel_length));
        vessel_segment_array[2+3*units_in_x_direction+2]->GetNode(1)->SetLocation(ChastePoint<DIM>(2*diagonal_vessel_length+horizontal_vessel_length, 2*diagonal_vessel_length));

        for (unsigned i = (2+3*units_in_x_direction+3); i < (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i++)
        {
                vessel_segment_array[i]->GetNode(1)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i-3]->GetNode(1)->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
                                vessel_segment_array[i-3]->GetNode(1)->GetLocation()[1]));
        }

        for (unsigned i = (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i < number_of_vessels - units_in_x_direction; i++)
        {
                vessel_segment_array[i]->GetNode(1)->SetLocation(ChastePoint<DIM>(vessel_segment_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode(1)->GetLocation()[0],
                                vessel_segment_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode(1)->GetLocation()[1] + 2*diagonal_vessel_length));
        }

        vessel_segment_array[number_of_vessels - units_in_x_direction]->GetNode(1)->SetLocation(ChastePoint<DIM>(2*(diagonal_vessel_length+horizontal_vessel_length),
                        vessel_segment_array[number_of_vessels - units_in_x_direction-1]->GetNode(0)->GetLocation()[1] + diagonal_vessel_length));

        for (unsigned i = 1; i < units_in_x_direction; i++)
        {
                vessel_segment_array[number_of_vessels - units_in_x_direction + i]->GetNode(1)->SetLocation(ChastePoint<DIM>(
                                vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(1)->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
                                vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(1)->GetLocation()[1]));
        }

        typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iterator;

        for (unsigned i = 0; i < vessel_segment_array.size(); i++)
        {
            vessel_array.push_back(CaVessel<DIM>::Create(vessel_segment_array[i]));
        }

        pVesselNetwork->AddVessels(vessel_array);

        pVesselNetwork->MergeCoincidentNodes();

        return pVesselNetwork;
}


template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateNetworkFromVtkFile(const std::string& filename)
{
	// Create an empty vessel network
	boost::shared_ptr<CaVascularNetwork<DIM>  > pVesselNetwork(new CaVascularNetwork<DIM>());

	// Create a VTK PolyData object based on the contents of the input VTK file
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkPolyData* pPolyData = reader->GetOutput();

	std::vector<double> radii;

	// Create a vector of chaste points corresponding to vessel node locations
	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
	for(vtkIdType i = 0; i < pPolyData->GetNumberOfPoints(); i++)
	{
		double point_coords[3];
		pPolyData->GetPoint(i, point_coords);
		if (DIM < 3)
		{
			nodes.push_back(boost::shared_ptr<VascularNode<DIM> > (VascularNode<DIM>::Create(ChastePoint<DIM>(point_coords[0], point_coords[1]))));
		}
		else
		{
			nodes.push_back(boost::shared_ptr<VascularNode<DIM> > (VascularNode<DIM>::Create(ChastePoint<DIM>(point_coords[0], point_coords[1], point_coords[2]))));
		}
	}

	///\ todo Radii also can be input in the form of cell data, the reader should be able to
	// also handle this format. It would be nicer to just read everything into a map of string
	// to vector doubles as then other data could be read in without much modification.

	// Extract radii corresponding to each node from the VTK Polydata and store them in a list.
	vtkPointData* pPointData = pPolyData->GetPointData();
	std::string radius_label = "Radius";
	for(vtkIdType i = 0; i < pPointData->GetNumberOfArrays(); i++)
	{
		std::string array_name = pPointData->GetArrayName(i);
		if(array_name.compare(radius_label) == 0)
		{
            vtkDoubleArray* pScalars = vtkDoubleArray::SafeDownCast(pPointData->GetArray(i));
			for(vtkIdType j = 0; j < pScalars->GetNumberOfTuples(); j++)
			{
				radii.push_back(pScalars->GetValue(j));
			}
		}
	}

	// Extract vessels from the VTK Polydata. This is done by iterating over a VTK CellArray object which
	// returns a 'pointList' vtkIdList object. This object contains the point IDs of the nodes which make up
	// the vessel.
	vtkCellArray* pCellArray = pPolyData->GetLines();

	///\ todo Deleting this pointer at the end of the block caused problems, not sure why. Should check best way to delete
	// vtk pointers. Smart pointers don't seem to work for vtkIdType in v5.2.
	vtkIdType* pSegmentList;
	vtkIdType num_segments;

	for(int i = 0; i < pPolyData->GetNumberOfLines(); i++)
	{
		// Make a new vessel
		pCellArray->GetNextCell(num_segments, pSegmentList);

		///\ todo For now get the average radius over the whole vessel, would be good to get radii in individual
		// segments, otherwise we are loosing data from the vtk file.
		double average_radius = radii[pSegmentList[0]];

		std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
		// Add segments to the vessels in order
		for (int j = 1; j < num_segments; j++)
		{
			// make a new segment
			segments.push_back(boost::shared_ptr<CaVesselSegment<DIM> > (CaVesselSegment<DIM>::Create(nodes[pSegmentList[j-1]], nodes[pSegmentList[j]])));
			average_radius = average_radius + radii[pSegmentList[j]];
		}
		boost::shared_ptr<CaVessel<DIM> > vessel(CaVessel<DIM>::Create(segments));
		vessel->GetDataContainer().SetData("radius", average_radius/double(num_segments));

		// Add the resulting vessel to the network
		pVesselNetwork->AddVessel(vessel);
	}

	return pVesselNetwork;
}

//Explicit instantiation
template class VasculatureGenerator<2>;
template class VasculatureGenerator<3>;
