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

#ifndef SEGMENTFLOWPROPERTIES_HPP_
#define SEGMENTFLOWPROPERTIES_HPP_

#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>

/**
 * This is a class for vascular segment flow properties.
 *
 * This class stores segment data for vessel network flow problems. Each segment has
 * an instance of the class.
 */
class SegmentFlowProperties : public boost::enable_shared_from_this<SegmentFlowProperties>
{

private:

    /**
     * Haematocrit in the vessel at this segment
     */
    double mHaematocrit;

    /**
     * Blood flow rate in the vessel at this segment
     */
    double mFlowRate;

    /**
     * Impedance of this vessel segment
     */
    double mImpedance;

    /**
     * Viscosity of this vessel segment
     */
    double mViscosity;

    /**
     * Wall shear stress of this vessel segment
     */
    double mWallShearStress;

    /**
     * Mechanical stimulus of this vessel segment
     */
    double mMechanicalStimulus;

    /**
     * Metabolic stimulus of this vessel segment
     */
    double mMetabolicStimulus;

    /**
     * Upstream conducted stimulus of this vessel segment
     */
    double mUpstreamConductedStimulus;

    /**
     * Downstream conducted stimulus of this vessel segment
     */
    double mDownstreamConductedStimulus;

    /**
     * Shrinking stimulus of this vessel segment
     */
    double mShrinkingStimulus;

public:

    /**
     * Constructor
     *
     * @param pressure the pressure in the vessel at the node
     * @param isInputNode whether the node is an input node
     * @param isOutputNode whether the node is an output node
     */
    SegmentFlowProperties();

    /**
     * Destructor
     */
    ~SegmentFlowProperties();

    /**
     * Return the haematocrit
     *
     * @return the segment haematocrit
     */
    double GetHaematocrit() const;

    /**
     * Return the impedance
     *
     * @return the segment impedance
     */
    double GetImpedance() const;

    /**
     * Return the flow rate
     *
     * @return the segment flow rate
     */
    double GetFlowRate() const;

    /**
     * Return the segment flow velocity
     *
     * @return the segment velocity
     */
    double GetFlowVelocity() const;

    /**
     * Return the segment viscosity
     *
     * @return the segment viscosity
     */
    double GetViscosity() const;

    /**
     * Return the segment wall shear stress of this vessel segment
     *
     * @return the segment wall shear stress
     */
    double GetWallShearStress() const;

    /**
     * Return the segment mechanical stimulus of this vessel segment
     *
     * @return the segment mechanical stimulus
     */
    double GetMechanicalStimulus() const;

    /**
     * Return the metabolic stimulus of this vessel segment
     *
     * @return the segment metabolic stimulus
     */
    double GetMetabolicStimulus() const;

    /**
     * Return the upstream conducted stimulus of this vessel segment
     *
     * @return the segment upstream conducted stimulus
     */
    double GetUpstreamConductedStimulus() const;

    /**
     * Return the downstream conducted stimulus of this vessel segment
     *
     * @return the segment downstream conducted stimulus
     */
    double GetDownstreamConductedStimulus() const;

    /**
     * Return the shrinking stimulus of this vessel segment
     *
     * @return the segment shrinking stimulus
     */
    double GetShrinkingStimulus() const;

    /**
     * Return a map of segment data for use by the vtk writer
     *
     * @return a map of segment data for use by the vtk writer
     */
    std::map<std::string, double> GetVtkData() const;

    /**
     * Set the haematocrit
     *
     * @param haematocrit the haematocrit in the segment
     */
    void SetHaematocrit(double haematocrit);

    /**
     * Set the flow rate
     *
     * @param flowRate the flow rate in the segment
     */
    void SetFlowRate(double flowRate);

    /**
     * Set the impedance
     *
     * @param impedance the impedance in the segment
     */
    void SetImpedance(double impedance);

    /**
     * Set the viscosity
     *
     * @param viscosity the viscosity in the segment
     */
    void SetViscosity(double viscosity);

    /**
     * Set the wall shear stress of this vessel segment
     *
     * @param wallShear the wall shear stress in the segment
     */
    void SetWallShearStress(double wallShear);

    /**
     * Set the mechanical stimulus of this vessel segment
     *
     * @param mechStimulus the mechanical stimulus in the segment
     */
    void SetMechanicalStimulus(double mechStimulus);

    /**
     * Set the metabolic stimulus of this vessel segment
     *
     * @param metabolicStimulus the metabolic stimulus in the segment
     */
    void SetMetabolicStimulus(double metabolicStimulus);

    /**
     * Set the upstream conducted stimulus of this vessel segment
     *
     * @param upstreamStimulus the upstream conducted stimulus in the segment
     */
    void SetUpstreamConductedStimulus(double upstreamStimulus);

    /**
     * Set the downstream conducted stimulus of this vessel segment
     *
     * @param downstreamStimulus the downstream conducted stimulus in the segment
     */
    void SetDownstreamConductedStimulus(double downstreamStimulus);

    /**
     * Set the shrinking stimulus of this vessel segment
     *
     * @param shrinkingStimulus the shrinking stimulus in the segment
     */
    void SetShrinkingStimulus(double shrinkingStimulus);

};

#endif /* SEGMENTFLOWPROPERTIES_HPP_ */
