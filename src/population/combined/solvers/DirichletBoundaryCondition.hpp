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
#include "AbstractCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "RegularGrid.hpp"
#include "TetrahedralMesh.hpp"

/**
 * Helper struct for defining the type of boundary condition.
 */
struct BoundaryConditionType
{
    enum Value
    {
        POINT, FACET, OUTER, VESSEL_LINE, VESSEL_VOLUME, CELL, IN_PART
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

    /**
     * The vessel network, used for VESSEL_ type conditions
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /**
     * The cell population, used for CELL type conditions
     */
    boost::shared_ptr<AbstractCellPopulation<DIM> > mpCellPopulation;

    /**
     * A part for prescribing part and facet based conditions
     */
    boost::shared_ptr<Part<DIM> > mpDomain;

    /**
     * Point locations for POINT type conditions
     */
    std::vector<c_vector<double, DIM> > mPoints;

    /**
     * The type of boundary condition
     */
    BoundaryConditionType::Value mType;

    /**
     * Where the boundary condition value is obtained from
     */
    BoundaryConditionSource::Value mSource;

    /**
     * A label specifying the array name from which to obtain the condition magnitude. Used for LABEL
     * conditions.
     */
    std::string mLabel;

    /**
     * The prescribed value of the boundary condition.
     */
    double mValue;

    /**
     * The grid for solvers using regular grids
     */
    boost::shared_ptr<RegularGrid<DIM, DIM> > mpRegularGrid;

    /**
     * The mesh for solvers using finite element meshes
     */
    boost::shared_ptr<TetrahedralMesh<DIM, DIM> > mpMesh;

public:

    /**
     * Constructor
     */
    DirichletBoundaryCondition();

    /**
     * Destructor
     */
    virtual ~DirichletBoundaryCondition();

    /**
     * Factory constructor method
     */
    static boost::shared_ptr<DirichletBoundaryCondition<DIM> > Create();

    BoundaryConditionType::Value GetType();

    double GetValue();

    std::pair<bool, double> GetValue(c_vector<double,DIM> location, double tolerance);

    void UpdateBoundaryConditionContainer(boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> > pContainer);

    void UpdateRegularGridBoundaryConditions(boost::shared_ptr<std::vector<std::pair<bool, double> > > pBoundaryConditions, double tolerance);

    /**
     * Set the cell population, used in CELL type sources
     * @param rCellPopulation a reference to the cell population
     */
    void SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation);

    void SetDomain(boost::shared_ptr<Part<DIM> > pDomain);

    void SetLabelName(const std::string& label);

    /**
     * Set the points for POINT type sources
     * @param points the point locations for POINT type sources
     */
    void SetPoints(std::vector<c_vector<double, DIM> > points);

    /**
     * Set the finite element mesh
     * @param pMesh the finite element mesh
     */
    void SetMesh(boost::shared_ptr<TetrahedralMesh<DIM, DIM> > pMesh);

    /**
     * Set the regular grid
     * @param pRegularGrid the regular grid
     */
    void SetRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pRegularGrid);

    void SetSource(BoundaryConditionSource::Value boundarySource);

    void SetType(BoundaryConditionType::Value boundaryType);

    void SetValue(double value);

    /**
     * Set the vessel network used in VESSEL type sources
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);
};

#endif /* DIRICHLETBOUNDARYCONDITION_HPP_ */
