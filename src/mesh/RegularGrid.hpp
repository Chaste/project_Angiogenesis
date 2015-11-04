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
#include "UblasVectorInclude.hpp"
#include "SmartPointers.hpp"

/**
 * A simple description of a regular lattice for use in hybrid simulations.
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

public:

    /**
     *  Constructor
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
     * Get the 1-D grid index for given x,y,z indices
     * @param x_index the grid x index
     * @param y_index the grid y index
     * @param z_index the grid z index
     * @return the grid 1-d index
     */
    unsigned Get1dGridIndex(unsigned x_index, unsigned y_index, unsigned z_index);

    /* Return the grid extents in x, y, z
     * @return the grid extents
     */
    std::vector<unsigned> GetExtents();

    /*
     * Get the location of a point on the grid for given x,y,z indices
     */
    c_vector<double, SPACE_DIM> GetLocation(unsigned x_index, unsigned y_index, unsigned z_index);

    /*
     * Get the location of a point on the grid for given 1-d grid index
     */
    c_vector<double, SPACE_DIM> GetLocationOf1dIndex(unsigned grid_index);

    /* Return the origin in x, y, z
     * @return the grid origin
     */
    c_vector<double, SPACE_DIM> GetOrigin();

    unsigned GetNumberOfPoints();

    /**
     * Return the grid spacing
     * @return the grid spacing
     */
    double GetSpacing();

    bool IsOnBoundary(unsigned grid_index);

    bool IsOnBoundary(unsigned x_index, unsigned y_index, unsigned z_index);

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
};

#endif /* REGULARGRID_HPP_*/
