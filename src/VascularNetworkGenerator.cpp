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

#include "VascularNetworkGenerator.hpp"

template<unsigned DIM>
VascularNetworkGenerator<DIM>::VascularNetworkGenerator(boost::shared_ptr<CaVessel<DIM> > prototype)
	: mpPrototypeVessel(prototype),
	  mArterialHaematocritLevel(0.0),
	  mVenousOutputPressure(0.0),
	  mArterialInputPressure(0.0),
	  mInitialRadius(0.0)
{
}

template<unsigned DIM>
VascularNetworkGenerator<DIM>::~VascularNetworkGenerator()
{
}

template<unsigned DIM>
void VascularNetworkGenerator<DIM>::SetArterialHaematocritLevel(double value)
{
	mArterialHaematocritLevel = value;

	if(mArterialHaematocritLevel > 1.0 || mArterialHaematocritLevel < 0.0)
	{
		WARN_ONCE_ONLY("Specified Haematocrit is outside the physical range of 0 to 1, this may cause problems in subsequent calculations.");
	}
}

template<unsigned DIM>
void VascularNetworkGenerator<DIM>::SetArterialInputPressure(double value)
{
	mArterialInputPressure = value;

	if(mArterialInputPressure > 1.01*pow(10.0,5) || mArterialInputPressure < 1.01*pow(10.0,5)/76.0)
	{
		WARN_ONCE_ONLY("Specified Arterial Input Pressure is outside the physical range of 13.29 to 101.0e3, this may cause problems in subsequent calculations.");
	}
}

template<unsigned DIM>
void VascularNetworkGenerator<DIM>::SetVenousOutputPressure(double value)
{
	mVenousOutputPressure = value;

	if(mVenousOutputPressure > 1.01*pow(10.0,5) || mVenousOutputPressure < 1.01*pow(10.0,5)/76.0)
	{
		WARN_ONCE_ONLY("Specified Venous Output Pressure is outside the physical range of 13.29 to 101.0e3, this may cause problems in subsequent calculations.");
	}
}

template<unsigned DIM>
void VascularNetworkGenerator<DIM>::SetInitialRadius(double value)
{
	mInitialRadius = value;

	if(mInitialRadius > 50.0*pow(10.0,-6) || mInitialRadius < 1.0*pow(10.0,-6))
	{
		////todo 50 micron seems like a low upper limit to vessel radius, can this be increased?
		WARN_ONCE_ONLY("Specified Initial Radius is outside the physical range of 1 to 50 micron, this may cause problems in subsequent calculations.");
	}
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VascularNetworkGenerator<DIM>::GenerateHexagonalNetwork(unsigned width,
		unsigned height,
		unsigned vessel_length,
		std::string venous_output_position)
{

	// Ensure a suitable venous output position is specified
	if(venous_output_position != "North East" && venous_output_position != "South East")
	{
		EXCEPTION("Unrecognized venous output position specified, options are North East or SouthEast");
	}

	// Determine the number of repeating hexagon units in the x and y directions based on the input height, width and vessel length
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
	for (unsigned i = 0; i < number_of_vessels; i++)
	{
		vessel_array.push_back(CreateVessel());
	}

	// Set the left side coordinates of vessels in the network
	vessel_array[0]->SetNode1Location(ChastePoint<DIM>(0.0, 0.0));
	vessel_array[1]->SetNode1Location(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length,diagonal_vessel_length));
	vessel_array[2]->SetNode1Location(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length, 0.0));
	for (unsigned i = 3; i < (2+3*units_in_x_direction); i++)
	{
		vessel_array[i]->SetNode1Location(ChastePoint<DIM>(vessel_array[i-3]->GetNode1()->GetLocation()[0] +
				2.0*(diagonal_vessel_length+horizontal_vessel_length), vessel_array[i-3]->GetNode1()->GetLocation()[1]));
	}
	vessel_array[2+3*units_in_x_direction]->SetNode1Location(ChastePoint<DIM>(0, 2.0*diagonal_vessel_length));
	vessel_array[2+3*units_in_x_direction+1]->SetNode1Location(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
	vessel_array[2+3*units_in_x_direction+2]->SetNode1Location(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length, diagonal_vessel_length));

	for (unsigned i = (2+3*units_in_x_direction+3); i < (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i++)
	{
		vessel_array[i]->SetNode1Location(ChastePoint<DIM>(vessel_array[i-3]->GetNode1()->GetLocation()[0] +
				2.0*(diagonal_vessel_length+horizontal_vessel_length), vessel_array[i-3]->GetNode1()->GetLocation()[1]));
	}

	for (unsigned i = (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i < number_of_vessels - units_in_x_direction; i++)
	{
		vessel_array[i]->SetNode1Location(ChastePoint<DIM>(vessel_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode1()->GetLocation()[0],
				vessel_array[i-(2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode1()->GetLocation()[1] + 2.0*diagonal_vessel_length));
	}
	vessel_array[number_of_vessels-units_in_x_direction]->SetNode1Location(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length,
			vessel_array[number_of_vessels-units_in_x_direction-1]->GetNode1()->GetLocation()[1]+diagonal_vessel_length));

	for (unsigned i = 1; i < units_in_x_direction; i++)
	{
		vessel_array[number_of_vessels-units_in_x_direction+i]->SetNode1Location(ChastePoint<DIM>(
				vessel_array[number_of_vessels-units_in_x_direction+i-1]->GetNode1()->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
				vessel_array[number_of_vessels-units_in_x_direction+i-1]->GetNode1()->GetLocation()[1]));
	}

	// Set the right side coordinates of vessels in the network
	vessel_array[0]->SetNode2Location(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
	vessel_array[1]->SetNode2Location(ChastePoint<DIM>(2.0*diagonal_vessel_length+horizontal_vessel_length, 0.0));
	vessel_array[2]->SetNode2Location(ChastePoint<DIM>(2.0*(diagonal_vessel_length+horizontal_vessel_length), 0.0));

	for (unsigned i = 3; i < (2+3*units_in_x_direction); i++)
	{
		vessel_array[i]->SetNode2Location(ChastePoint<DIM>(vessel_array[i-3]->GetNode2()->GetLocation()[0]+2*(diagonal_vessel_length+horizontal_vessel_length),
				vessel_array[i-3]->GetNode2()->GetLocation()[1]));
	}
	vessel_array[2+3*units_in_x_direction]->SetNode2Location(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
	vessel_array[2+3*units_in_x_direction+1]->SetNode2Location(ChastePoint<DIM>(diagonal_vessel_length+horizontal_vessel_length, diagonal_vessel_length));
	vessel_array[2+3*units_in_x_direction+2]->SetNode2Location(ChastePoint<DIM>(2*diagonal_vessel_length+horizontal_vessel_length, 2*diagonal_vessel_length));

	for (unsigned i = (2+3*units_in_x_direction+3); i < (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i++)
	{
		vessel_array[i]->SetNode2Location(ChastePoint<DIM>(vessel_array[i-3]->GetNode2()->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
				vessel_array[i-3]->GetNode2()->GetLocation()[1]));
	}

	for (unsigned i = (2+3*units_in_x_direction+3*(units_in_x_direction+1)); i < number_of_vessels - units_in_x_direction; i++)
	{
		vessel_array[i]->SetNode2Location(ChastePoint<DIM>(vessel_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode2()->GetLocation()[0],
				vessel_array[i - (2+3*units_in_x_direction+3*(units_in_x_direction+1))]->GetNode2()->GetLocation()[1] + 2*diagonal_vessel_length));
	}

	vessel_array[number_of_vessels - units_in_x_direction]->SetNode2Location(ChastePoint<DIM>(2*(diagonal_vessel_length+horizontal_vessel_length),
			vessel_array[number_of_vessels - units_in_x_direction-1]->GetNode1()->GetLocation()[1] + diagonal_vessel_length));

	for (unsigned i = 1; i < units_in_x_direction; i++)
	{
		vessel_array[number_of_vessels - units_in_x_direction + i]->SetNode2Location(ChastePoint<DIM>(
				vessel_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode2()->GetLocation()[0] + 2*(diagonal_vessel_length+horizontal_vessel_length),
				vessel_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode2()->GetLocation()[1]));
	}

	// Create  a list of vessel segment coordinates for each vessel in network
	for(unsigned i = 0; i < number_of_vessels; i++)
	{
		if (vessel_array[i]->GetNode2()->GetLocation()[1] == vessel_array[i]->GetNode1()->GetLocation()[1])
		{
			// Vessel is horizontal
			for(unsigned j = 0; j < horizontal_vessel_length + 1; j++)
			{
				vessel_array[i]->SetNextVesselSegmentCoordinate(ChastePoint<DIM>(vessel_array[i]->GetNode1()->GetLocation()[0]+j,
						vessel_array[i]->GetNode1()->GetLocation()[1]));
			}
		}
		else if (vessel_array[i]->GetNode2()->GetLocation()[1] > vessel_array[i]->GetNode1()->GetLocation()[1])
		{
			// Vessel is diagonal - going downwards
			for(unsigned j = 0; j < diagonal_vessel_length + 1; j++)
			{
				vessel_array[i]->SetNextVesselSegmentCoordinate(ChastePoint<DIM>(vessel_array[i]->GetNode1()->GetLocation()[0]+j,
						vessel_array[i]->GetNode1()->GetLocation()[1]+j));
			}
		}
		else if (vessel_array[i]->GetNode2()->GetLocation()[1] < vessel_array[i]->GetNode1()->GetLocation()[1])
		{
			// Vessel is diagonal - going upwards
			for(unsigned j = 0; j < diagonal_vessel_length + 1; j++)
			{
				vessel_array[i]->SetNextVesselSegmentCoordinate(ChastePoint<DIM>(vessel_array[i]->GetNode1()->GetLocation()[0]+j,
						vessel_array[i]->GetNode1()->GetLocation()[1]-j));
			}
		}
	}

	// Add the vessels to the network
	for (unsigned i=0; i < number_of_vessels; i++)
	{
		pVesselNetwork->AddVessel(vessel_array[i]);
	}

	// Initialize vessel lengths, radii and haematocrit
	for (unsigned i = 0; i < pVesselNetwork->GetNumberOfVesselsInNetwork(); i++)
	{
		pVesselNetwork->GetVessel(i)->SetRadius(mInitialRadius);
		pVesselNetwork->GetVessel(i)->SetHaematocritLevel(mArterialHaematocritLevel);
	}
	pVesselNetwork->SetArterialHaematocritLevel(mArterialHaematocritLevel);
	pVesselNetwork->SetArterialInputPressure(mArterialInputPressure);
	pVesselNetwork->SetVenousOutputPressure(mVenousOutputPressure);

	if (venous_output_position == "North East")
	{
		pVesselNetwork->SetOutputNode(pVesselNetwork->GetVessel(2+3*units_in_x_direction-1)->GetNode2()->GetLocation());
	}
	else if (venous_output_position == "South East")
	{
		pVesselNetwork->SetOutputNode(pVesselNetwork->GetVessel(pVesselNetwork->GetNumberOfVesselsInNetwork()-1-units_in_x_direction)->GetNode2()->GetLocation());
	}

	// Set the input node
	pVesselNetwork->SetInputNode(ChastePoint<DIM>(0.0, 0.0));

	return pVesselNetwork;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VascularNetworkGenerator<DIM>::GenerateSingleBifurcationNetwork(unsigned vessel_length)
{

	// Make the vessels



}

#ifdef CHASTE_VTK
template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VascularNetworkGenerator<DIM>::GenerateNetworkFromVtkFile(std::string filename)
{
	// Create an empty vessel network
	boost::shared_ptr<CaVascularNetwork<DIM>  > pVesselNetwork(new CaVascularNetwork<DIM>());

	// Create a VTK PolyData object based on the contents of the input VTK file
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkPolyData* pPolyData = reader->GetOutput();
	std::vector<ChastePoint<DIM> > point_locations;
	std::vector<double> radii;

	// Create a vector of chaste points corresponding to vessel node locations
	for(vtkIdType i = 0; i < pPolyData->GetNumberOfPoints(); i++)
	{
		double point_coords[3];
		pPolyData->GetPoint(i, point_coords);
		if (DIM < 3)
		{
			point_locations.push_back(ChastePoint<DIM>(point_coords[0], point_coords[1]));
		}
		else
		{
			point_locations.push_back(ChastePoint<DIM>(point_coords[0], point_coords[1], point_coords[2]));
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
		boost::shared_ptr<CaVessel<DIM> > vessel = CreateVessel();
		pCellArray->GetNextCell(num_segments, pSegmentList);
		vessel->SetNode1Location(point_locations[pSegmentList[0]]);
		vessel->SetNode2Location(point_locations[pSegmentList[num_segments - 1]]);

		///\ todo For now get the average radius over the whole vessel, would be good to get radii in individual
		// segments, otherwise we are loosing data from the vtk file.
		double average_radius = radii[pSegmentList[0]];

		// Add segments to the vessels in order
		for (int j = 0; j < num_segments; j++)
		{
			vessel->SetNextVesselSegmentCoordinate(point_locations[pSegmentList[j]]);
			average_radius = average_radius + radii[pSegmentList[j]];
		}
		vessel->SetRadius(average_radius/double(num_segments));
		vessel->SetHaematocritLevel(mArterialHaematocritLevel);

		// Add the resulting vessel to the network
		pVesselNetwork->AddVessel(vessel);
	}

	pVesselNetwork->SetArterialHaematocritLevel(mArterialHaematocritLevel);
	pVesselNetwork->SetArterialInputPressure(mArterialInputPressure);
	pVesselNetwork->SetVenousOutputPressure(mVenousOutputPressure);

	return pVesselNetwork;
}
#endif // CHASTE_VTK

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > VascularNetworkGenerator<DIM>::CreateVessel()
{
	boost::shared_ptr<CaVessel<DIM> > vessel(new CaVessel<DIM>());

	// If there is a prototype vessel use it to make sure that all vessels in the network "know about" the right chemicals and
	// that those chemical species have the right concentrations.
	// Check if pointer has been assigned
	if(mpPrototypeVessel)
	{
		for (unsigned i = 0; i < mpPrototypeVessel->GetNumberOfIntraVascularChemicals(); i++)
		{
			// \todo This is aweful! Need to sort out the interfaces for intra vascular chemicals and collections to make them easier to handle.
			vessel->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical(mpPrototypeVessel->GetCollectionOfIntraVascularChemicals().
					GetIntraVascularChemicalCollection()[i].GetChemicalName(),
					Concentration(mpPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetConcentration(),
							mpPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetUnits()),
							mpPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetPermeabilityOfVesselWallToChemical());
		}
	}
	// we wish to count vessels instantiated now as part of existing vasculature so that they can be
	// distinguished from vessels forming the neovasculature when those vessels form
	vessel->SetIsPartOfNeovasculature(false);

	return vessel;
}

//Explicit instantiation
template class VascularNetworkGenerator<2>;
template class VascularNetworkGenerator<3>;
