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
#include "UnitCollection.hpp"
#include "VesselSegment.hpp"
#include "SmartPointers.hpp"
#include "AbstractVesselNetworkComponentFlowProperties.hpp"
#include "SimulationTime.hpp"

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
     * Return the dimensionless haematocrit
     * @return the segment haematocrit
     */
    double GetHaematocrit(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the 'dimensional' haematocrit
     * @return the segment haematocrit
     */
    units::quantity<unit::dimensionless> GetDimensionalHaematocrit(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the haematocrit in SI units
     * @return the segment haematocrit
     */
    double GetHaematocritSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the impedance
     *
     * @return the segment impedance
     */
    double GetImpedance(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the impedance
     *
     * @return the segment impedance
     */
    units::quantity<unit::flow_impedance> GetDimensionalImpedance(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the impedance
     *
     * @return the segment impedance
     */
    double GetImpedanceSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the flow rate
     *
     * @return the segment flow rate
     */
    double GetFlowRate(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the flow rate
     *
     * @return the segment flow rate
     */
    units::quantity<unit::flow_rate> GetDimensionalFlowRate(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the flow rate
     *
     * @return the segment flow rate
     */
    double GetFlowRateSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment flow velocity
     *
     * @return the segment velocity
     */
    double GetFlowVelocity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment flow velocity
     *
     * @return the segment velocity
     */
    units::quantity<unit::velocity> GetDimensionalFlowVelocity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment flow velocity
     *
     * @return the segment velocity
     */
    double GetFlowVelocitySI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment viscosity
     *
     * @return the segment viscosity
     */
    double GetViscosity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment viscosity
     *
     * @return the segment viscosity
     */
    units::quantity<unit::dynamic_viscosity> GetDimensionalViscosity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment viscosity
     *
     * @return the segment viscosity
     */
    double GetViscositySI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment wall shear stress of this vessel segment
     *
     * @return the segment wall shear stress
     */
    double GetWallShearStress(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment wall shear stress of this vessel segment
     *
     * @return the segment wall shear stress
     */
    units::quantity<unit::pressure> GetDimensionalWallShearStress(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the segment wall shear stress of this vessel segment
     *
     * @return the segment wall shear stress
     */
    double GetWallShearStressSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the growth stimulus of this vessel segment
     *
     * @return the segment growth stimulus
     */
    double GetGrowthStimulus(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the growth stimulus of this vessel segment
     *
     * @return the segment growth stimulus
     */
    units::quantity<unit::rate> GetDimensionalGrowthStimulus(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

    /**
     * Return the growth stimulus of this vessel segment
     *
     * @return the segment growth stimulus
     */
    double GetGrowthStimulusSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const;

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
    bool HasVesselRegressed();

    /**
     * Rescue vessel from regression.
     */
    void ResetRegressionTimer();

    /**
     * Set the haematocrit
     *
     * @param haematocrit the haematocrit in the segment
     */
    void SetDimensionalHaematocrit(units::quantity<unit::dimensionless> haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the haematocrit
     *
     * @param haematocrit the haematocrit in the segment
     */
    void SetHaematocrit(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the haematocrit
     *
     * @param haematocrit the haematocrit in the segment
     */
    void SetHaematocritSI(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the flow rate
     *
     * @param flowRate the flow rate in the segment
     */
    void SetFlowRate(double flowRate, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the flow rate
     *
     * @param flowRate the flow rate in the segment
     */
    void SetDimensionalFlowRate(units::quantity<unit::flow_rate> flowRate, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the flow rate
     *
     * @param flowRate the flow rate in the segment
     */
    void SetFlowRateSI(double flowRate, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the impedance
     *
     * @param impedance the impedance in the segment
     */
    void SetImpedance(double impedance, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the impedance
     *
     * @param impedance the impedance in the segment
     */
    void SetDimensionalImpedance(units::quantity<unit::flow_impedance> impedance, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the impedance
     *
     * @param impedance the impedance in the segment
     */
    void SetImpedanceSI(double impedance, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the viscosity
     *
     * @param viscosity the viscosity in the segment
     */
    void SetDimensionalViscosity(units::quantity<unit::dynamic_viscosity> viscosity, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the viscosity
     *
     * @param viscosity the viscosity in the segment
     */
    void SetViscosity(double viscosity, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the viscosity
     *
     * @param viscosity the viscosity in the segment
     */
    void SetViscositySI(double viscosity, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the wall shear stress of this vessel segment
     *
     * @param wallShear the wall shear stress in the segment
     */
    void SetDimensionalWallShearStress(units::quantity<unit::pressure> wallShear, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the wall shear stress of this vessel segment
     *
     * @param wallShear the wall shear stress in the segment
     */
    void SetWallShearStress(double wallShear, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the wall shear stress of this vessel segment
     *
     * @param wallShear the wall shear stress in the segment
     */
    void SetWallShearStressSI(double wallShear, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the growth stimulus of this vessel segment
     *
     * @param mechStimulus the growth stimulus in the segment
     */
    void SetGrowthStimulus(double stimulus, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the growth stimulus of this vessel segment
     *
     * @param mechStimulus the growth stimulus in the segment
     */
    void SetGrowthStimulusSI(double stimulus, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set the growth stimulus of this vessel segment
     *
     * @param mechStimulus the growth stimulus in the segment
     */
    void SetDimensionalGrowthStimulus(units::quantity<unit::rate> stimulus, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments);

    /**
     * Set dimensionless time until removal of vessel from network.
     * @param time the time until vessel removal
     */
    void SetTimeUntilRegression(double time);

    /**
     * Set dimensional time until removal of vessel from network.
     * @param time the time until vessel removal
     */
    void SetDimensionalTimeUntilRegression(units::quantity<unit::time> time);

    /**
     * Set time until removal of vessel from network.
     * @param time the time until vessel removal
     */
    void SetTimeUntilRegressionSI(double time);

};

#endif /* VesselFlowProperties_HPP_ */
