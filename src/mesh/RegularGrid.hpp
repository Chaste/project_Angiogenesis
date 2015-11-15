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

#ifndef REGULARGRID_HPP_
#define REGULARGRID_HPP_

#include <vector>
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellPopulation.hpp"
#include "Part.hpp"

/**
 * A class for describing regular grids, calculating point and line to grid point relationships and
 * storing cell and vessel to grid point maps.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class RegularGrid
{
    /**
     *  The spacing between grid points
     */
    double mSpacing;

    /**
     *  The number of grid points in each direction
     */
    std::vector<unsigned> mExtents;

    /**
     *  The origin of the grid in x,y,z. Corresponds to location of front, bottom, left corner.
     */
    c_vector<double, SPACE_DIM> mOrigin;

    /**
     * The vessel network
     */
    boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > mpNetwork;

    /**
     * The cell population. This memory pointed to is not managed in this class.
     */
    AbstractCellPopulation<SPACE_DIM>* mpCellPopulation;

    /**
     * A map of cells corresponding to a point on the grid
     */
    std::vector<std::vector<CellPtr> > mPointCellMap;

    /**
     * A map of vessel segments corresponding to a point on the grid
     */
    std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > > mPointSegmentMap;

public:

    /**
     * Constructor
     */
    RegularGrid();

    /**
     * Factory constructor method
     * @return a shared pointer to a new grid
     */
    static boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > Create();

    /**
     * Desctructor
     */
    ~RegularGrid();

    /**
     * Generate a grid based on the bounding box of the supplied part
     * @param pPart the part from which to get the bounding box
     * @param gridSize the grid spacing
     */
    void GenerateFromPart(boost::shared_ptr<Part<SPACE_DIM> > pPart, double gridSize);

    /**
     * Get the 1-D grid index for given x,y,z indices
     * @param x_index the grid x index
     * @param y_index the grid y index
     * @param z_index the grid z index
     * @return the grid 1-d index
     */
    unsigned Get1dGridIndex(unsigned xIndex, unsigned yIndex, unsigned zIndex);

    /**
     * Return the grid extents in x, y, z. Always dimension 3.
     * @return the grid extents
     */
    std::vector<unsigned> GetExtents();

    /**
     * Get the location of a point on the grid for given x, y ,z indices
     * @param x_index the grid x index
     * @param y_index the grid y index
     * @param z_index the grid z index
     * @return the location of the point
     */
    c_vector<double, SPACE_DIM> GetLocation(unsigned xIndex, unsigned yIndex, unsigned zIndex);

    /*
     * Get the location of a point on the grid for given 1-d grid index
     * @param gridIndex the 1d grid index
     * @return the location of the point
     */
    c_vector<double, SPACE_DIM> GetLocationOf1dIndex(unsigned gridIndex);

    /**
     * Get all of the grid locations
     * @return a vector containing all grid locations in grid order
     */
    std::vector<c_vector<double, SPACE_DIM> > GetLocations();

    /**
     * Return the number of points in the grid
     * @return the number of points in the grid
     */
    unsigned GetNumberOfPoints();

    /**
     * Return the origin in x, y, z
     * @return the grid origin
     */
    c_vector<double, SPACE_DIM> GetOrigin();

    /**
     * Return a vector of input point indices which in the bounding boxes of each grid point
     * @bool inputPoints a vector of point locations
     * @return the indices of input points in the bounding box of each grid point
     */
    std::vector<std::vector<unsigned> > GetPointPointMap(std::vector<c_vector<double, SPACE_DIM> > inputPoints);

    /**
     * Return the point cell map
     * @bool update update the map
     * @return the point cell map
     */
    const std::vector<std::vector<CellPtr> >& GetPointCellMap(bool update = true);

    /**
     * Return the point segments map
     * @bool update update the map
     * @return the point segment map
     */
    std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > > GetPointSegmentMap(bool update = true, bool useVesselSurface = false);

    /**
     * Return the grid spacing
     * @return the grid spacing
     */
    double GetSpacing();

    std::vector<double> InterpolateGridValues(std::vector<c_vector<double, SPACE_DIM> > locations, std::vector<double> values, bool useVtk = false);


	//double VolumeAverageQuantity(std::vector<double> values);
    /**
     * Is the input location in the bounding box of the grid point
     * @param point the location of interest
     * @param gridIndex the grid point of interest
     * @return is the input location in the bounding box of the grid point
     */
    bool IsLocationInPointVolume(c_vector<double, SPACE_DIM> point, unsigned gridIndex);

    /**
     * Is the point on the outer boundary of the domain
     * @param the 1d grid index gridIndex
     * @return is the point on the outer boundary of the domain
     */
    bool IsOnBoundary(unsigned gridIndex);

    /**
     * Is the point on the outer boundary of the domain
     * @param xIndex the grid x index
     * @param yIndex the grid y index
     * @param zIndex the grid z index
     * @return is the point on the outer boundary of the domain
     */
    bool IsOnBoundary(unsigned xIndex, unsigned yIndex, unsigned zIndex);

    /**
     * Set the cell population
     * @param rCellPopulation a reference to the cell population
     */
    void SetCellPopulation(AbstractCellPopulation<SPACE_DIM>& rCellPopulation);

    /* Set the grid extents in x, y, z
     * @param extents the grid extents
     */
    void SetExtents(std::vector<unsigned> extents);

    /* Set the origin in x, y, z
     * @param origin the grid origin
     */
    void SetOrigin(c_vector<double, SPACE_DIM> origin);

    /**
     * Set the grid spacing
     * @param spacing the grid spacing
     */
    void SetSpacing(double spacing);

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > pNetwork);
};

#endif /* REGULARGRID_HPP_*/
