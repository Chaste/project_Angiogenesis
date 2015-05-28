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

#ifndef NODEFLOWPROPERTIES_HPP_
#define NODEFLOWPROPERTIES_HPP_

#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>

/**
 * This is a class for vascular node flow properties.
 *
 * This class stores nodal data for vessel network flow problems. Each node has
 * an instance of the class.
 */
class NodeFlowProperties : public boost::enable_shared_from_this<NodeFlowProperties>
{

private:

    /**
     * Pressure in the vessel at this node
     */
    double mPressure;

    /**
     *  Is the node an input node
     */
    bool mIsInputNode;

    /**
     *  Is the node an output node
     */
    bool mIsOutputNode;


public:

    /**
     * Constructor
     */
    NodeFlowProperties();

    /**
     * Destructor
     */
    ~NodeFlowProperties();

    /**
     * Return the pressure in the vessel at the node
     *
     * @return the pressure in the vessel at the node
     */
    double GetPressure() const;

    /**
     * Return a map of nodal data for use by the vtk writer
     *
     * @return a map of nodal data for use by the vtk writer
     */
    std::map<std::string, double> GetVtkData() const;

    /**
     * Return true if the node is an input node
     *
     * @return true if the node is an input node
     */
    bool IsInputNode() const;

    /**
     * Return true if the node is an output node
     *
     * @return true if the node is an output node
     */
    bool IsOutputNode() const;

    /**
     * Set the vessel pressure at this node
     *
     * @param pressure the vessel pressure
     */
    void SetPressure(double pressure);

    /**
     * Set that the node is an input node
     * @param isInput whether the node is an input
     */
    void SetIsInputNode(bool isInput);

    /**
     *  Set that the node is an output node
     * @param isOutput whether the node is an output
     */
    void SetIsOutputNode(bool isOutput);
};

#endif /* NODEFLOWPROPERTIES_HPP_ */
