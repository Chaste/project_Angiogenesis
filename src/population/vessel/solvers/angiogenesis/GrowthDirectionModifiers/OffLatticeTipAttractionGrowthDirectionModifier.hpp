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

#ifndef OFFLATTICETDGROWTHDIRECTIONMODIFIER_HPP_
#define OFFLATTICETDGROWTHDIRECTIONMODIFIER_HPP_

#include <vector>
#include <string>
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "AbstractGrowthDirectionModifier.hpp"
#include "CaVascularNetwork.hpp"
#include "VascularNode.hpp"
#include "UblasIncludes.hpp"

template<unsigned DIM>
class OffLatticeTipAttractionGrowthDirectionModifier : public AbstractGrowthDirectionModifier<DIM>
{

private:

    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

public:

    /**
     * Constructor.
     */
    OffLatticeTipAttractionGrowthDirectionModifier();

    /**
     * Destructor.
     */
    virtual ~OffLatticeTipAttractionGrowthDirectionModifier();

    static boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<DIM> > Create();

    void SetNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);
    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */


    c_vector<double, DIM> GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                             boost::shared_ptr<VascularNode<DIM> > pNode);

};

#endif /* OFFLATTICETDGROWTHDIRECTIONMODIFIER_HPP_ */
