/*
 * VascularNetworkGenerator.cpp
 *
 *  Created on: 8 Jan 2015
 *      Author: chaste
 */

#include "VascularNetworkGenerator.hpp"

VascularNetworkGenerator::VascularNetworkGenerator()
{

}

VascularNetworkGenerator::~VascularNetworkGenerator()
{

}

std::vector<unsigned> VascularNetworkGenerator::GenerateHexagonallyTesselatedVascularNetwork(PottsMesh<2>* p_mesh, unsigned width, unsigned height)
{

    int horizontalVesselLength = 5;
    int diagonalVesselLength = 5;
    int numberOfRepeatingUnitsInXDirection = floor((((width) - horizontalVesselLength - 2*diagonalVesselLength - 1)/(2*(horizontalVesselLength+diagonalVesselLength))));
    int numberOfRepeatingUnitsInYDirection = floor((((height)-1)/(2*diagonalVesselLength)));

}
