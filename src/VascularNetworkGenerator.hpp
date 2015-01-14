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
     * @param p_mesh Potts mesh
     * @param width width of Potts mesh
     * @param height width of Potts mesh
     * @return a vector of location indices which specify the locations of vascular cells inside a
     *         hexagonally tesselated vessel network
     */
    boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > GenerateHexagonallyTesselatedVascularNetwork(PottsMesh<2>* p_mesh, unsigned width, unsigned height);

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
