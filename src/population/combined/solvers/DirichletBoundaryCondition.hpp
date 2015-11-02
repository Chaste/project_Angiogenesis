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

#ifndef DIRICHLETBOUNDARYCONDITION_HPP_
#define DIRICHLETBOUNDARYCONDITION_HPP_

#include <vector>
#include <string>
#include "UblasIncludes.hpp"
//#include "SimpleCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"

/*
 * Helper struct for defining the type of boundary condition.
 * It can be, point, multipoint, facet, vessel-line, vessel-volume, cell-point or cell-volume
 */
struct BoundaryConditionType
{
    enum Value
    {
        POINT, MULTI_POINT, OUTER, OUTER_2D, FACET, VESSEL_LINE, VESSEL_VOLUME, CELL_POINT, IN_PART
    };
};

/*
 * Helper struct for defining the source of the boundary condition value.
 * It can be from a data array or a single prescribed value.
 */
struct BoundaryConditionSource
{
    enum Value
    {
        LABEL_BASED, PRESCRIBED
    };
};

/*
 * An class for describing dirichlet boundary conditions for use with some hybrid solvers.
 */

template<unsigned DIM>
class DirichletBoundaryCondition
{
private:

    /* The vessel network
    */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /* The cell population
    */
//    boost::shared_ptr<SimpleCellPopulation<DIM> > mpCellPopulation;

    boost::shared_ptr<Part<DIM> > mpDomain;

    std::vector<c_vector<double, DIM> > mPoints;

    BoundaryConditionType::Value mType;

    BoundaryConditionSource::Value mSource;

    std::string mLabel;

    double mValue;

public:

    /* Constructor
     */
    DirichletBoundaryCondition();

    /* Destructor
     */
    virtual ~DirichletBoundaryCondition();

    /* Factory constructor method
     */
    static boost::shared_ptr<DirichletBoundaryCondition<DIM> > Create();

    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

//    void SetCellPopulation(boost::shared_ptr<SimpleCellPopulation<DIM> > pCellPopulation);

    void SetDomain(boost::shared_ptr<Part<DIM> > pDomain);

    void SetPoint(c_vector<double, DIM> point);

    void SetPoints(std::vector<c_vector<double, DIM> > points);

    void SetType(BoundaryConditionType::Value boundaryType);

    BoundaryConditionType::Value GetType();

    void SetSource(BoundaryConditionSource::Value boundarySource);

    void SetLabelName(const std::string& label);

    void SetValue(double value);

    std::pair<bool, double> GetValue(c_vector<double, DIM> location = zero_vector<double>(DIM), double tolerance = 1.e-6);

    std::vector<std::pair<bool, double> > GetValues(std::vector<c_vector<double, DIM> > locations, double tolerance = 1.e-6);
};

#endif /* DIRICHLETBOUNDARYCONDITION_HPP_ */
