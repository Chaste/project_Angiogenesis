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

#ifndef GREENSFUNCTIONSOLVER_HPP_
#define GREENSFUNCTIONSOLVER_HPP_

#include <vector>
#include <string>
#include <map>
#include <boost/multi_array.hpp>
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "CaVesselSegment.hpp"
#include "Part.hpp"
#include "AbstractRegularGridHybridSolver.hpp"

template<unsigned DIM>
class GreensFunctionSolver : public AbstractRegularGridHybridSolver<DIM>
{
    boost::shared_ptr<Part<DIM> > mpDomain;
    std::vector<ChastePoint<DIM> > mSinkCoordinates;
    std::vector<unsigned> mSinkPointMap;
    std::vector<ChastePoint<DIM> > mSubSegmentCoordinates;
    std::vector<double> mSubSegmentLengths;
    std::vector<double> mSinkRates;
    std::vector<double> mSourceRates;
    std::vector<double> mSegmentConcentration;
    std::vector<double> mTissueConcentration;
    std::map<unsigned, boost::shared_ptr<CaVesselSegment<DIM> > > mSegmentPointMap;
    boost::shared_ptr<boost::multi_array<double, 2> > mGtt;
    boost::shared_ptr<boost::multi_array<double, 2> > mGvv;
    boost::shared_ptr<boost::multi_array<double, 2> > mGvt;
    boost::shared_ptr<boost::multi_array<double, 2> > mGtv;

public:

    GreensFunctionSolver();

    ~GreensFunctionSolver();

    void Solve(bool writeSolution = false);

private:

    void GenerateSubSegments();

    void GenerateTissuePoints();

    void UpdateGreensFunctionMatrices(bool updateGtt = 0, bool updateGvv = 0, bool updateGtv = 0, bool updateGvt = 0);

    boost::shared_ptr<boost::multi_array<double, 2> > GetVesselVesselInteractionMatrix();

    boost::shared_ptr<boost::multi_array<double, 2> > GetTissueTissueInteractionMatrix();

    boost::shared_ptr<boost::multi_array<double, 2> > GetTissueVesselInteractionMatrix();

    boost::shared_ptr<boost::multi_array<double, 2> > GetVesselTissueInteractionMatrix();

    void WriteSolution(std::map<std::string, std::vector<double> >& segmentPointData);
};

#endif /* GREENSFUNCTIONSOLVER_HPP_ */
