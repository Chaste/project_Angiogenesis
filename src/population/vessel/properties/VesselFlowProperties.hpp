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

#ifndef VesselFlowProperties_HPP_
#define VesselFlowProperties_HPP_

#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>

#include "BaseUnits.hpp"
#include "UnitCollection.hpp"
#include "VesselSegment.hpp"
#include "SmartPointers.hpp"
#include "AbstractVesselNetworkComponentFlowProperties.hpp"

/**
 * This is a class for vessel flow properties.
 *
 * This class calculates vessel data for vessel network flow problems. Each vessel has
 * an instance of the class.
 */
template<unsigned DIM>
class VesselFlowProperties : public boost::enable_shared_from_this<VesselFlowProperties<DIM> >, public AbstractVesselNetworkComponentFlowProperties<DIM>
{

    using AbstractVesselNetworkComponentFlowProperties<DIM>::GetOutputData;

private:

    /**
     * Whether a vessel is currently undergoing regression. A vessel can be saved from this fate.
     */
    bool mUndergoingRegression;

    /**
     * Whether a vessel should be removed from the network. A vessel exists inside the network until they are removed.
     */
    bool mRemoveViaRegression;

    /**
     * When the vessel will be removed
     */
    units::quantity<unit::time> mRegressionTime;

public:

    /**
     * Constructor
     *
     * @param pressure the pressure in the vessel at the node
     * @param isInputNode whether the node is an input node
     * @param isOutputNode whether the node is an output node
     */
    VesselFlowProperties();

    /**
     * Destructor
     */
    ~VesselFlowProperties();

    /**
     * Return the 'dimensional' haematocrit
     * @return the segment haematocrit
     */
    units::quantity<unit::dimensionless> GetHaematocrit(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the impedance
     *
     * @return the segment impedance
     */
    units::quantity<unit::flow_impedance> GetImpedance(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the flow rate
     *
     * @return the segment flow rate
     */
    units::quantity<unit::flow_rate> GetFlowRate(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment viscosity
     *
     * @return the segment viscosity
     */
    units::quantity<unit::dynamic_viscosity> GetViscosity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment wall shear stress of this vessel segment
     *
     * @return the segment wall shear stress
     */
    units::quantity<unit::pressure> GetWallShearStress(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the growth stimulus of this vessel segment
     *
     * @return the segment growth stimulus
     */
    units::quantity<unit::rate> GetGrowthStimulus(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return a map of segment data for use by the vtk writer
     *
     * @return a map of segment data for use by the vtk writer
     */
    std::map<std::string, double> GetOutputData(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * @return whether regression timer has started.
     */
    bool HasRegressionTimerStarted();

    /**
     * @return whether the vessel should regress (be removed).
     */
    bool HasVesselRegressed(units::quantity<unit::time> simulationReferenceTime);

    /**
     * Rescue vessel from regression.
     */
    void ResetRegressionTimer();

    /**
     * Set the haematocrit
     *
     * @param haematocrit the haematocrit in the segment
     */
    void SetHaematocrit(units::quantity<unit::dimensionless> haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the flow rate
     *
     * @param flowRate the flow rate in the segment
     */
    void SetFlowRate(units::quantity<unit::flow_rate> flowRate, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the impedance
     *
     * @param impedance the impedance in the segment
     */
    void SetImpedance(units::quantity<unit::flow_impedance> impedance, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the viscosity
     *
     * @param viscosity the viscosity in the segment
     */
    void SetViscosity(units::quantity<unit::dynamic_viscosity> viscosity, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the wall shear stress of this vessel segment
     *
     * @param wallShear the wall shear stress in the segment
     */
    void SetWallShearStress(units::quantity<unit::pressure> wallShear, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the growth stimulus of this vessel segment
     *
     * @param mechStimulus the growth stimulus in the segment
     */
    void SetGrowthStimulus(units::quantity<unit::rate> stimulus, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set dimensional time until removal of vessel from network.
     * @param time the time until vessel removal
     */
    void SetTimeUntilRegression(units::quantity<unit::time> time, units::quantity<unit::time> simulationReferenceTime);

};

#endif /* VesselFlowProperties_HPP_ */
