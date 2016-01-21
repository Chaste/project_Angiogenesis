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
#include <map>
#include "UblasIncludes.hpp"
#include "RegularGrid.hpp"
#include "HybridMesh.hpp"

/**
 * Specify the type of source.
 * POINT: Source locations are given by a vector of c_vectors
 * VESSEL: Source locations are along a vessel
 * CELL: Source locations are cell centres
 * SOLUTION: Source strength depends on a previous solution at this location
 */
struct SourceType
{
    enum Value
    {
        POINT, VESSEL, CELL, SOLUTION
    };
};

/**
 * Specify whether a single prescribed value is used for the source (PRESCRIBED) or
 * whether the source strength depends on a specified label or map (LABEL).
 */
struct SourceStrength
{
    enum Value
    {
        LABEL, PRESCRIBED
    };
};

/**
 * This class calculates the value of discrete sources at grid/mesh locations in continuum problems.
 * A grid or mesh is set and the source values are calculated at each grid point or in each element.
 * Child classes can be used to customize the way source strengths are calculated.
 */
template<unsigned DIM>
class DiscreteSource
{

protected:

    /**
     * The grid for solvers using regular grids
     */
    boost::shared_ptr<RegularGrid<DIM, DIM> > mpRegularGrid;

    /**
     * The mesh for the finite element solver
     */
    boost::shared_ptr<HybridMesh<DIM, DIM> > mpMesh;

    /**
     * A continuum field sampled on the regular grid points or mesh nodes
     */
    std::vector<double>  mpSolution;

    /**
     * Locations for POINT type sources
     */
    std::vector<c_vector<double, DIM> > mPoints;

    /**
     * The type of source
     */
    SourceType::Value mType;

    /**
     * Where the source strength is obtained from
     */
    SourceStrength::Value mSourceStrength;

    /**
     * A label specifying the array name from which to obtain the source strength. Used for LABEL
     * source strengths.
     */
    std::string mLabel;

    /**
     * The prescribed value of the source strength. Used for PRESCRIBED source strengths.
     */
    double mValue;

    /**
     * Is the source linear in the solution, if not it is constant in it.
     */
    bool mIsLinearInSolution;

public:

    /**
     *  Constructor
     */
    DiscreteSource();

    /**
     * Destructor
     */
    virtual ~DiscreteSource();

    /**
     * Factory constructor method
     * @return a pointer to an instance of the class
     */
    static boost::shared_ptr<DiscreteSource<DIM> > Create();

    /**
     * Return the type of source, (POINT, VESSEL, etc.)
     * @return an enum giving the type of source
     */
    SourceType::Value GetType();

    /**
     * Return the values of the source strengths sampled on the mesh elements
     * @return a vector of source strengths
     */
    std::vector<double> GetMeshValues();

    /**
     * Return the values of the source strengths sampled on the regular grid
     * @return a vector of source strengths
     */
    std::vector<double> GetRegularGridValues();

    /**
     * Return the values of the POINT source strengths sampled on the regular grid
     * @return a vector of source strengths
     */
    virtual std::vector<double> GetPointRegularGridValues();

    /**
     * Return the values of the CELL source strengths sampled on the regular grid
     * @return a vector of source strengths
     */
    virtual std::vector<double> GetCellRegularGridValues();

    /**
     * Return the values of the VESSEL source strengths sampled on the regular grid
     * @return a vector of source strengths
     */
    virtual std::vector<double> GetVesselRegularGridValues();

    /**
     * Return the values of the SOLUTION source strengths sampled on the regular grid
     * @return a vector of source strengths
     */
    virtual std::vector<double> GetSolutionDependentRegularGridValues();

    /**
     * Is the source strength linear in the solution variable, otherwise it is constant
     * @return a bool asking is the source strength linear in the solution variable
     */
    bool IsLinearInSolution();

    /**
     * Set whether the source strength is linear in the solution
     * @param isLinear a bool asking is the source strength linear in the solution variable
     */
    void SetIsLinearInSolution(bool isLinear);

    /**
     * Set the name of the label used in LABEL type sources
     * @param rLabel the label for the source strength value
     */
    void SetLabelName(const std::string& rLabel);

    /**
     * Set the points for POINT type sources
     * @param points the point locations for POINT type sources
     */
    void SetPoints(std::vector<c_vector<double, DIM> > points);

    /**
     * Set the regular grid
     * @param pRegularGrid the regular grid
     */
    void SetRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pRegularGrid);

    /**
     * Set the finite element mesh
     * @param pMesh the finite element mesh
     */
    void SetMesh(boost::shared_ptr<HybridMesh<DIM, DIM> > pMesh);

    /**
     * Set the sampled field from which to obtain a solution for SOLUTION type sources
     * @param pSolution the field from which to use solution values
     */
    void SetSolution(std::vector<double> solution);

    /**
     * Set where the value of the source strength is obtained, e.g. LABEL, PRESCRIBED
     * @param boundarySource enum specifying where the value of the source strength is obtained
     */
    void SetSource(SourceStrength::Value boundarySource);

    /**
     * Set the type of source, e.g. CELL, VESSEL
     * @param boundaryType enum specifying the type of source
     */
    void SetType(SourceType::Value boundaryType);

    /**
     * Set the value of the source for PRESCRIBED type sources
     * @param value the value of the source
     */
    void SetValue(double value);
};

#endif /* DISCRETESOURCE_HPP_ */
