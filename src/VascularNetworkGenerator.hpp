/*
 * VascularNetworkGenerator.hpp
 *
 *  Created on: 8 Jan 2015
 *      Author: chaste
 */

#ifndef VASCULARNETWORKGENERATOR_HPP_
#define VASCULARNETWORKGENERATOR_HPP_

class VascularNetworkGenerator
{

public:

    /*
     * Constructor
     */
    VascularNetworkGenerator();

    /*
     * Destructor
     */
    ~VascularNetworkGenerator();

    /*
     * Returns the location indices for vascular cells in a hexagonally tesselated vessel network
     *
     * @param p_mesh Potts mesh
     * @param width width of Potts mesh
     * @param height width of Potts mesh
     * @return a vector of location indices which specify the locations of vascular cells inside a
     *         hexagonally tesselated vessel network
     */
    std::vector<unsigned> GenerateHexagonallyTesselatedVascularNetwork(PottsMesh<2>* p_mesh, unsigned width, unsigned height);

};

#endif /* VASCULARNETWORKGENERATOR_HPP_ */
