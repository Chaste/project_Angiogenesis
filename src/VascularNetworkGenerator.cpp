/*
 * VascularNetworkGenerator.cpp
 *
 *  Created on: 8 Jan 2015
 *      Author: chaste
 */

#include "VascularNetworkGenerator.hpp"
#include <string>
#include <vector>

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
boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > VascularNetworkGenerator<SPATIAL_DIM>::GenerateHexagonallyTesselatedVascularNetwork(PottsMesh<2>* p_mesh, unsigned width, unsigned height)
{

    std::string VenousOutputPosition = "North East";

    double horizontalVesselLength = 5;
    double diagonalVesselLength = 5;
    int numberOfRepeatingUnitsInXDirection = floor((((width) - horizontalVesselLength - 2*diagonalVesselLength - 1)/(2*(horizontalVesselLength+diagonalVesselLength))));
    int numberOfRepeatingUnitsInYDirection = floor((((height)-1)/(2*diagonalVesselLength)));

    int numberOfVessels = (numberOfRepeatingUnitsInXDirection*numberOfRepeatingUnitsInYDirection*6) + numberOfRepeatingUnitsInYDirection*5 + numberOfRepeatingUnitsInXDirection;
    int numberOfNodes = (4*(numberOfRepeatingUnitsInXDirection+1))*(numberOfRepeatingUnitsInYDirection) + 2*(numberOfRepeatingUnitsInXDirection+1);


    boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > VN(new CaVascularNetwork<SPATIAL_DIM>());

    std::vector<boost::shared_ptr<CaVessel<SPATIAL_DIM> > > vesselArray(numberOfVessels);

    for (int i = 0; i < numberOfVessels; i++)
    {
        vesselArray[i] = CreateVessel();
    }

    // initialize left side coordinates of vessels in network

    vesselArray[0]->SetNode1Location(ChastePoint<SPATIAL_DIM>(0,0));
    vesselArray[1]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));
    vesselArray[2]->SetNode1Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,0));

    for (int i = 3; i < (2+3*numberOfRepeatingUnitsInXDirection); i++)
    {
        vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode1()->GetLocation()[1]));
    }

    vesselArray[2+3*numberOfRepeatingUnitsInXDirection]->SetNode1Location(ChastePoint<SPATIAL_DIM>(0,2*diagonalVesselLength));
    vesselArray[2+3*numberOfRepeatingUnitsInXDirection+1]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
    vesselArray[2+3*numberOfRepeatingUnitsInXDirection+2]->SetNode1Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));

    for (int i = (2+3*numberOfRepeatingUnitsInXDirection+3); i < (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i++)
    {
        vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode1()->GetLocation()[1]));
    }

    for (int i = (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i < numberOfVessels - numberOfRepeatingUnitsInXDirection; i++)
    {

        vesselArray[i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode1()->GetLocation()[0],vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode1()->GetLocation()[1] + 2*diagonalVesselLength));
    }

    vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection]->SetNode1Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection - 1]->GetNode1()->GetLocation()[1] + diagonalVesselLength));


    for (int i = 1; i < numberOfRepeatingUnitsInXDirection; i++)
    {
        vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i]->SetNode1Location(ChastePoint<SPATIAL_DIM>(vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode1()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode1()->GetLocation()[1]));
    }

    // initialize Right side coordinates of vessels in network

    vesselArray[0]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
    vesselArray[1]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,0));
    vesselArray[2]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*(diagonalVesselLength+horizontalVesselLength),0));


    for (int i = 3; i < (2+3*numberOfRepeatingUnitsInXDirection); i++)
    {
        vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode2()->GetLocation()[1]));
    }
    vesselArray[2+3*numberOfRepeatingUnitsInXDirection]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength,diagonalVesselLength));
    vesselArray[2+3*numberOfRepeatingUnitsInXDirection+1]->SetNode2Location(ChastePoint<SPATIAL_DIM>(diagonalVesselLength+horizontalVesselLength,diagonalVesselLength));
    vesselArray[2+3*numberOfRepeatingUnitsInXDirection+2]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*diagonalVesselLength+horizontalVesselLength,2*diagonalVesselLength));

    for (int i = (2+3*numberOfRepeatingUnitsInXDirection+3); i < (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i++)
    {
        vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i-3]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[i-3]->GetNode2()->GetLocation()[1]));
    }

    for (int i = (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1)); i < numberOfVessels - numberOfRepeatingUnitsInXDirection; i++)
    {
        vesselArray[i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode2()->GetLocation()[0],vesselArray[i - (2+3*numberOfRepeatingUnitsInXDirection+3*(numberOfRepeatingUnitsInXDirection+1))]->GetNode2()->GetLocation()[1] + 2*diagonalVesselLength));
    }


    vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection]->SetNode2Location(ChastePoint<SPATIAL_DIM>(2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection-1]->GetNode1()->GetLocation()[1] + diagonalVesselLength));


    for (int i = 1; i < numberOfRepeatingUnitsInXDirection; i++)
    {

        vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i]->SetNode2Location(ChastePoint<SPATIAL_DIM>(vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode2()->GetLocation()[0] + 2*(diagonalVesselLength+horizontalVesselLength),vesselArray[numberOfVessels - numberOfRepeatingUnitsInXDirection + i - 1]->GetNode2()->GetLocation()[1]));
    }

    // create list of vessel segment coordinates for each vessel in network

    for( int i = 0; i < numberOfVessels; i++)
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


    for (int i=0; i < numberOfVessels; i++)
    {
        VN->AddVessel(vesselArray[i]);
    }

    // initialize vessel lengths, radii and haematocrit

    for (int i = 0; i < VN->GetNumberOfVesselsInNetwork(); i++)
    {
        VN->GetVessel(i)->SetRadius(mInitialRadius);
        VN->GetVessel(i)->SetHaematocritLevel(mArterialHaematocritLevel);
    }

    VN->SetArterialHaematocritLevel(mArterialHaematocritLevel);
    VN->SetArterialInputPressure(mArterialInputPressure);
    VN->SetVenousOutputPressure(mVenousOutputPressure);



    return VN;
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > VascularNetworkGenerator<SPATIAL_DIM>::CreateVessel()
{

    boost::shared_ptr<CaVessel<SPATIAL_DIM> > V(new CaVessel<SPATIAL_DIM>());

    // use prototype vessel to make sure that all vessels in the network "know about" the right chemicals and
    // that those chemical species have the right concentrations.

    for (int i = 0; i < pPrototypeVessel->GetNumberOfIntraVascularChemicals(); i++)
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
