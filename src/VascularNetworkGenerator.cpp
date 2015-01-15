/*
 * VascularNetworkGenerator.cpp
 *
 *  Created on: 8 Jan 2015
 *      Author: chaste
 */

#include "VascularNetworkGenerator.hpp"
#include <string>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>

template<unsigned SPATIAL_DIM>
VascularNetworkGenerator<SPATIAL_DIM>::VascularNetworkGenerator(boost::shared_ptr<CaVessel<SPATIAL_DIM> > prototype) : pPrototypeVessel(prototype), mArterialHaematocritLevel(0.45), mVenousOutputPressure(15*(1.01*pow(10.0,5)/760)), mArterialInputPressure(25*(1.01*pow(10.0,5)/760)), mInitialRadius(10*pow(10.0,-6))
{

}

template<unsigned SPATIAL_DIM>
VascularNetworkGenerator<SPATIAL_DIM>::~VascularNetworkGenerator()
{

}

template<unsigned SPATIAL_DIM>
void VascularNetworkGenerator<SPATIAL_DIM>::SetArterialHaematocritLevel(double value)
{
	mArterialHaematocritLevel = value;
	assert(mArterialHaematocritLevel < 1);
	assert(mArterialHaematocritLevel > 0);
}

template<unsigned SPATIAL_DIM>
void VascularNetworkGenerator<SPATIAL_DIM>::SetArterialInputPressure(double value)
{
	mArterialInputPressure = value;
	assert(mArterialInputPressure > 0.1*(1.01*pow(10.0,5)/760));
	assert(mArterialInputPressure < 760*(1.01*pow(10.0,5)/760));
}

template<unsigned SPATIAL_DIM>
void VascularNetworkGenerator<SPATIAL_DIM>::SetVenousOutputPressure(double value)
{
	mVenousOutputPressure = value;
	assert(mVenousOutputPressure > 0.1*(1.01*pow(10.0,5)/760));
	assert(mVenousOutputPressure < 760*(1.01*pow(10.0,5)/760));
}

template<unsigned SPATIAL_DIM>
void VascularNetworkGenerator<SPATIAL_DIM>::SetInitialRadius(double value)
{
	mInitialRadius = value;
	assert(mInitialRadius > 1*pow(10.0,-6));
	assert(mInitialRadius < 50*pow(10.0,-6));
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > VascularNetworkGenerator<SPATIAL_DIM>::GenerateHexagonallyTesselatedVascularNetwork(unsigned width,unsigned height,unsigned vessel_length,std::string venous_output_position)
{

	assert(venous_output_position == "North East" || venous_output_position == "South East");

	double horizontalVesselLength = vessel_length;
	double diagonalVesselLength = floor(vessel_length/std::sqrt(2));
	unsigned numberOfRepeatingUnitsInXDirection = floor((((width) - horizontalVesselLength - 2*diagonalVesselLength - 1)/(2*(horizontalVesselLength+diagonalVesselLength))));
	unsigned numberOfRepeatingUnitsInYDirection = floor((((height)-1)/(2*diagonalVesselLength)));

	assert(numberOfRepeatingUnitsInXDirection > 1);
	assert(numberOfRepeatingUnitsInYDirection > 1);


	unsigned numberOfVessels = (numberOfRepeatingUnitsInXDirection*numberOfRepeatingUnitsInYDirection*6) + numberOfRepeatingUnitsInYDirection*5 + numberOfRepeatingUnitsInXDirection;

	boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > VN(new CaVascularNetwork<SPATIAL_DIM>());

	std::vector<boost::shared_ptr<CaVessel<SPATIAL_DIM> > > vesselArray(numberOfVessels);

	for (unsigned i = 0; i < numberOfVessels; i++)
	{
		vesselArray[i] = CreateVessel();
	}

	// initialize left side coordinates of vessels in network

	vesselArray[0]->SetNode1Location(ChastePoint<SPATIAL_DIM>(0,0));
	vesselArray[1]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));
	vesselArray[2]->SetNode1Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,0));

	for (unsigned i = 3; i < (2+3*numberOfRepeatingUnitsInXDirection); i++)
	{
		vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode1()->GetLocation()[1]));
	}

	vesselArray[2+3*numberOfRepeatingUnitsInXDirection]->SetNode1Location(ChastePoint<SPATIAL_DIM>(0,2*diagonalVesselLength));
	vesselArray[2+3*numberOfRepeatingUnitsInXDirection+1]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
	vesselArray[2+3*numberOfRepeatingUnitsInXDirection+2]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));

	for (unsigned i = (2+3*numberOfRepeatingUnitsInXDirection+3); i < (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i++)
	{
		vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode1()->GetLocation()[1]));
	}

	for (unsigned i = (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i < numberOfVessels - numberOfRepeatingUnitsInXDirection; i++)
	{

		vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode1()->GetLocation()[0],vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode1()->GetLocation()[1] + 2*diagonalVesselLength));
	}

	vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection]->SetNode1Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection - 1]->GetNode1()->GetLocation()[1] + diagonalVesselLength));


	for (unsigned i = 1; i < numberOfRepeatingUnitsInXDirection; i++)
	{
		vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode1()->GetLocation()[1]));
	}

	// initialize Right side coordinates of vessels in network

	vesselArray[0]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
	vesselArray[1]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,0));
	vesselArray[2]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*(diagonalVesselLength+horizontalVesselLength),0));


	for (unsigned i = 3; i < (2+3*numberOfRepeatingUnitsInXDirection); i++)
	{
		vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode2()->GetLocation()[1]));
	}
	vesselArray[2+3*numberOfRepeatingUnitsInXDirection]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
	vesselArray[2+3*numberOfRepeatingUnitsInXDirection+1]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));
	vesselArray[2+3*numberOfRepeatingUnitsInXDirection+2]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,2*diagonalVesselLength));

	for (unsigned i = (2+3*numberOfRepeatingUnitsInXDirection+3); i < (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i++)
	{
		vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode2()->GetLocation()[1]));
	}

	for (unsigned i = (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i < numberOfVessels - numberOfRepeatingUnitsInXDirection; i++)
	{
		vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode2()->GetLocation()[0],vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode2()->GetLocation()[1] + 2*diagonalVesselLength));
	}


	vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection-1]->GetNode1()->GetLocation()[1] + diagonalVesselLength));


	for (unsigned i = 1; i < numberOfRepeatingUnitsInXDirection; i++)
	{

		vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode2()->GetLocation()[1]));
	}

	// create list of vessel segment coordinates for each vessel in network

	for( unsigned i = 0; i < numberOfVessels; i++)
	{

		if (vesselArray[i]->GetNode2()->GetLocation()[1] == vesselArray[i]->GetNode1()->GetLocation()[1])
		{
			// vessel is horizontal

			for(int j = 0; j < horizontalVesselLength + 1; j++)
			{
				vesselArray[i]->SetNextVesselSegmentCoordinate(ChastePoint<SPATIAL_DIM>(vesselArray[i]->GetNode1()->GetLocation()[0] + j,vesselArray[i]->GetNode1()->GetLocation()[1]));
			}
		}

		if (vesselArray[i]->GetNode2()->GetLocation()[1] > vesselArray[i]->GetNode1()->GetLocation()[1])
		{
			// vessel is diagonal - going downwards

			for(int j = 0; j < diagonalVesselLength + 1; j++)
			{
				vesselArray[i]->SetNextVesselSegmentCoordinate(ChastePoint<SPATIAL_DIM>(vesselArray[i]->GetNode1()->GetLocation()[0] + j,vesselArray[i]->GetNode1()->GetLocation()[1] + j));
			}
		}

		if (vesselArray[i]->GetNode2()->GetLocation()[1] < vesselArray[i]->GetNode1()->GetLocation()[1])
		{
			// vessel is diagonal - going upwards

			for(int j = 0; j < diagonalVesselLength + 1; j++)
			{
				vesselArray[i]->SetNextVesselSegmentCoordinate(ChastePoint<SPATIAL_DIM>(vesselArray[i]->GetNode1()->GetLocation()[0] + j,vesselArray[i]->GetNode1()->GetLocation()[1] - j));
			}
		}


	}


	for (unsigned i=0; i < numberOfVessels; i++)
	{
		VN->AddVessel(vesselArray[i]);
	}

	// initialize vessel lengths, radii and haematocrit

	for (unsigned i = 0; i < VN->GetNumberOfVesselsInNetwork(); i++)
	{
		VN->GetVessel(i)->SetRadius(mInitialRadius);
		VN->GetVessel(i)->SetHaematocritLevel(mArterialHaematocritLevel);
	}

	VN->SetArterialHaematocritLevel(mArterialHaematocritLevel);
	VN->SetArterialInputPressure(mArterialInputPressure);
	VN->SetVenousOutputPressure(mVenousOutputPressure);

	if (venous_output_position == "North East")
	{
		VN->SetOutputNode(VN->GetVessel(2+3*numberOfRepeatingUnitsInXDirection - 1)->GetNode2()->GetLocation());
	}

	if (venous_output_position == "South East")
	{
		VN->SetOutputNode(VN->GetVessel(VN->GetNumberOfVesselsInNetwork() - 1 - numberOfRepeatingUnitsInXDirection)->GetNode2()->GetLocation());
	}

	VN->SetInputNode(ChastePoint<SPATIAL_DIM>(0,0));

	return VN;
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > VascularNetworkGenerator<SPATIAL_DIM>::GenerateVascularNetworkFromVtkFile(std::string filename)
{

	// Create an empty vessel network
	boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM>  > VN(new CaVascularNetwork<SPATIAL_DIM>());

	// Create a VTK PolyData object based on the contents of the input VTK file
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkPolyData *p_polydata = reader->GetOutput();
	std::vector<ChastePoint<SPATIAL_DIM> > m_pointLocations;
	std::vector<double> m_radii;

	std::cout << "1" << std::endl;

	// Create a Vector of spatial location objects corresponding to node locations.
	for(vtkIdType i = 0; i < p_polydata->GetNumberOfPoints(); i++)
	{
		if (SPATIAL_DIM < 3)
		{
			double pointCoords[3];
			p_polydata->GetPoint(i,pointCoords);
			m_pointLocations.push_back(ChastePoint<SPATIAL_DIM>(pointCoords[0],pointCoords[1]));
		}
		else
		{
			double pointCoords[3];
			p_polydata->GetPoint(i,pointCoords);
			m_pointLocations.push_back(ChastePoint<SPATIAL_DIM>(pointCoords[0],pointCoords[1],pointCoords[2]));
		}
	}
	std::cout << "1" << std::endl;
	// Extract radii corresponding to each node from the VTK Polydata and store them in a list.
	vtkCellData *p_info = p_polydata->GetCellData();
	std::string radiusLabel = "Radius";
	for(vtkIdType i = 0; i < p_info->GetNumberOfArrays(); i++)
	{
	    std::cout << "number of arrays" << p_info->GetNumberOfArrays() << std::endl;
		std::string arrayName = p_info->GetArrayName(i);
		std::cout << arrayName << std::endl;
		if(arrayName.compare(radiusLabel) == 0)
		{
		    std::cout << "1" << std::endl;
            vtkDoubleArray* scalars = vtkDoubleArray::SafeDownCast(p_info->GetArray(i));
            std::cout << "2" << std::endl;

			for(vtkIdType j = 0; j < scalars->GetNumberOfTuples(); j++)
			{
			    std::cout << scalars->GetValue(j) << std::endl;

				m_radii.push_back(scalars->GetValue(j));
			}
		}
	}

	std::cout << "1" << std::endl;

	// Extract vessels from the VTK Polydata. This is done by iterating over a VTK CellArray object which
	// returns a 'pointList' vtkIdList object. This object contains the point IDs of the nodes which make up
	// the vessel.
	int numCells = p_polydata->GetNumberOfLines();
	vtkCellArray *p_cellArray = p_polydata->GetLines();
	vtkIdType* segList;
	vtkIdType num_segments;


	for(int i = 0; i < numCells; i++)
	{
		boost::shared_ptr<CaVessel<SPATIAL_DIM> > V = CreateVessel();

		p_cellArray->GetNextCell(num_segments, segList);

		// Set the vessel start and end nodes
		V->SetNode1Location(m_pointLocations[segList[0]]);
		V->SetNode2Location(m_pointLocations[segList[num_segments - 1]]);

		// Add segments to the vessels in order
		for (int j = 0; j < num_segments; j++)
		{
			V->SetNextVesselSegmentCoordinate(m_pointLocations[segList[j]]);
		}
		V->SetRadius(m_radii[i]);
		V->SetHaematocritLevel(mArterialHaematocritLevel);

		// Add the resulting vessel to the network
		VN->AddVessel(V);
	}

	VN->SetArterialHaematocritLevel(mArterialHaematocritLevel);
	VN->SetArterialInputPressure(mArterialInputPressure);
	VN->SetVenousOutputPressure(mVenousOutputPressure);

	delete  [] segList;

	return VN;

}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > VascularNetworkGenerator<SPATIAL_DIM>::CreateVessel()
{

	boost::shared_ptr<CaVessel<SPATIAL_DIM> > V(new CaVessel<SPATIAL_DIM>());

	// use prototype vessel to make sure that all vessels in the network "know about" the right chemicals and
	// that those chemical species have the right concentrations.

	for (unsigned i = 0; i < pPrototypeVessel->GetNumberOfIntraVascularChemicals(); i++)
	{
		// \todo This is aweful! Need to sort out the interfaces for intra vascular chemicals and collections to make them easier to handle.
		V->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical(pPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetChemicalName(),
				Concentration(pPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetConcentration(),
						pPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetUnits()),
						pPrototypeVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetPermeabilityOfVesselWallToChemical());
	}

	// we wish to count vessels instantiated now as part of existing vasculature so that they can be
	// distinguished from vessels forming the neovasculature when those vessels form
	V->SetIsPartOfNeovasculature(false);

	return V;

}



//Explicit instantiation

template class VascularNetworkGenerator<2>;
template class VascularNetworkGenerator<3>;
