/*
 * VascularNetworkGenerator.hpp
 *
 *  Created on: 8 Jan 2015
 *      Author: chaste
 */

#ifndef VASCULARNETWORKGENERATOR_HPP_
#define VASCULARNETWORKGENERATOR_HPP_

#include <vector>
#include "PottsMesh.hpp"
#include <boost/shared_ptr.hpp>

#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "CaVascularNetworkNode.hpp"

template<unsigned SPATIAL_DIM>
class VascularNetworkGenerator
{

private:

    /**
            Prototypes Vessel object from which other vessel objects will be instantiated.  This allows us to configure which chemicals
            a vessel must "know about" in each simulation.
     */
    boost::shared_ptr<CaVessel<SPATIAL_DIM> > pPrototypeVessel;

    /**
            The arterial haematocrit level that will be prescribed to the VesselNetwork.
     */
    double mArterialHaematocritLevel;

    /**
            The venous output pressure that will be prescribed to the VesselNetwork.
     */
    double mVenousOutputPressure;

    /**
            The arterial input pressure that will be prescribed to the VesselNetwork.
     */
    double mArterialInputPressure;

    /**
            The initial radius that will be prescribed to all Vessel objects in the VesselNewtork (unless radii
            are prescribed in a file, for example).
     */
    double mInitialRadius;

public:

    /*
     * Constructor
     *
     * @param prototype prototype vessel using which other vessel objects will be instantiated.
     */
    VascularNetworkGenerator(boost::shared_ptr<CaVessel<SPATIAL_DIM> > prototype);

    /*
     * Destructor
     */
    ~VascularNetworkGenerator();

    /*
     * Creates a hexagonally tesselated vessel network
     *
     * @param width width of Potts mesh.
     * @param height width of Potts mesh.
     * @param vessel_length approximate vessel length.
     * @param venous_output_position output node location on Potts mesh (can be "North East" or "South East").
     * @return a hexagonally tesselated vessel network
     */
    boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > GenerateHexagonallyTesselatedVascularNetwork(unsigned width,
    																				unsigned height,
    																				unsigned vessel_length,
    																				std::string venous_output_position);

    /*
         * Generates a vessel network from a vtk file.
         *
         * @param filename name of file in which vascular network is described.
         * @return the vascular network described in the prescribed file.
         */
        boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > GenerateVascularNetworkFromVtkFile(std::string filename);

    /**
       Assigns the arterial haematocrit level that will be prescribed to the VesselNetwork.

       @param value newly prescribed arterial haematocrit level.
     */
    void SetArterialHaematocritLevel(double value);

    /**
       Assigns the arterial input pressure that will be prescribed to the VesselNetwork.

       @param value newly prescribed arterial input pressure.
     */
    void SetArterialInputPressure(double value);

    /**
       Assigns the venous output pressure that will be prescribed to the VesselNetwork.

       @param value newly prescribed venous output pressure.
     */
    void SetVenousOutputPressure(double pressure);

    /**
       Assigns the initial radius that will be prescribed to all Vessel objects added to the VesselNetwork.

       @param value newly prescribed initial vessel radius.
     */
    void SetInitialRadius(double radius);

private:

    /*
     * Create a new instance of a basic prototype vessel.
     *
     * @return new instance of prototype vessel.
     */
    boost::shared_ptr<CaVessel<SPATIAL_DIM> >  CreateVessel();
};

#endif /* VASCULARNETWORKGENERATOR_HPP_ */
