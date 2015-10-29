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

#ifndef DISCRETESOURCE_HPP_
#define DISCRETESOURCE_HPP_

#include <vector>
#include <string>
#include <vector>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include "UblasIncludes.hpp"
#include "AbstractCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "AbstractCellProperty.hpp"
#include "ApoptoticCellProperty.hpp"

/*
 * Helper struct for defining the type of source.
 * It can be point, multipoint, facet, vessel-line, vessel-volume, cell-point or cell-volume
 */
struct SourceType
{
    enum Value
    {
        POINT, MULTI_POINT, VESSEL_LINE, CELL_POINT, SOLUTION
    };
};

/*
 * Helper struct for defining the source strength.
 * It can be from a data array or a single prescribed value.
 */
struct SourceStrength
{
    enum Value
    {
        LABEL_BASED, PRESCRIBED
    };
};


/*
 * An class for describing discrete sources for use with some hybrid solvers.
 */

template<unsigned DIM>
class DiscreteSource
{
private:

    /* The vessel network
    */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /* The cell population
    */
    boost::shared_ptr<AbstractCellPopulation<DIM> > mpCellPopulation;

    boost::shared_ptr<Part<DIM> > mpDomain;

    vtkSmartPointer<vtkImageData>  mpSolution;

    std::vector<c_vector<double, DIM> > mPoints;

    SourceType::Value mType;

    SourceStrength::Value mSourceStrength;

    std::string mLabel;

    double mValue;

    bool mIsLinearInSolution;

    std::vector<std::pair<AbstractCellProperty, double > > mMutationSpecificConsumptionRateMap;

public:

    /* Constructor
     */
    DiscreteSource();

    /* Destructor
     */
    virtual ~DiscreteSource();

    /* Factory constructor method
     */
    static boost::shared_ptr<DiscreteSource<DIM> > Create();

    void SetMutationSpecificConsumptionRateMap(std::vector<std::pair<AbstractCellProperty, double > > mutationSpecificConsumptionRateMap);

    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    void SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation);

    void SetDomain(boost::shared_ptr<Part<DIM> > pDomain);

    void SetPoint(c_vector<double, DIM> point);

    void SetPoints(std::vector<c_vector<double, DIM> > points);

    void SetSolution(vtkSmartPointer<vtkImageData>  pSolution);

    void SetType(SourceType::Value boundaryType);

    SourceType::Value GetType();

    void SetSource(SourceStrength::Value boundarySource);

    void SetLabelName(const std::string& label);

    bool IsLinearInSolution();

    void SetIsLinearInSolution(bool isLinear);

    void SetValue(double value);

    std::pair<bool, double> GetValue(c_vector<double, DIM> location = zero_vector<double>(DIM), double tolerance = 1.e-6);

    std::vector<std::pair<bool, double> > GetValues(std::vector<c_vector<double, DIM> > locations, double tolerance = 1.e-6);
};

#endif /* DISCRETESOURCE_HPP_ */
