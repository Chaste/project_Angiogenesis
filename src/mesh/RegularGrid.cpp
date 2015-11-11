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

#include "Exception.hpp"
#include "RegularGrid.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RegularGrid<ELEMENT_DIM, SPACE_DIM>::RegularGrid() :
    mSpacing(1.0),
    mExtents(std::vector<unsigned>(3,10)),
    mOrigin(zero_vector<double>(SPACE_DIM))
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::Create()
{
    typedef RegularGrid<ELEMENT_DIM, SPACE_DIM> Reg_Grid_Templated;
    MAKE_PTR(Reg_Grid_Templated, pSelf);
    return pSelf;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RegularGrid<ELEMENT_DIM, SPACE_DIM>::~RegularGrid()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RegularGrid<ELEMENT_DIM, SPACE_DIM>::Get1dGridIndex(unsigned x_index, unsigned y_index, unsigned z_index)
{
    return x_index + mExtents[0] * y_index + mExtents[0] * mExtents[1] * z_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetExtents()
{
    return mExtents;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocation(unsigned x_index, unsigned y_index, unsigned z_index)
{
    c_vector<double, SPACE_DIM> location;
    location[0] = double(x_index) * mSpacing + mOrigin[0];
    location[1] = double(y_index) * mSpacing + mOrigin[1];
    if(SPACE_DIM == 3)
    {
        location[2] = double(z_index) * mSpacing + mOrigin[2];
    }
    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocations()
{
    std::vector<c_vector<double, SPACE_DIM> > locations(GetNumberOfPoints());
    for(unsigned idx=0; idx<GetNumberOfPoints(); idx++)
    {
        locations[idx] = GetLocationOf1dIndex(idx);
    }
    return locations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool RegularGrid<ELEMENT_DIM, SPACE_DIM>::IsOnBoundary(unsigned grid_index)
{
    unsigned mod_z = grid_index % (mExtents[0] * mExtents[1]);
    unsigned z_index = (grid_index - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;
    return IsOnBoundary(x_index, y_index, z_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocationOf1dIndex(unsigned grid_index)
{
    unsigned mod_z = grid_index % (mExtents[0] * mExtents[1]);
    unsigned z_index = (grid_index - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;
    return GetLocation(x_index, y_index, z_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetOrigin()
{
    return mOrigin;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetNumberOfPoints()
{
    return mExtents[0] * mExtents[1] * mExtents[2];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetSpacing()
{
    return mSpacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool RegularGrid<ELEMENT_DIM, SPACE_DIM>::IsOnBoundary(unsigned x_index, unsigned y_index, unsigned z_index)
{
    if(x_index == 0 || x_index == mExtents[0] - 1)
    {
        return true;
    }
    if(y_index == 0 || y_index == mExtents[1] - 1)
    {
        return true;
    }
    if(ELEMENT_DIM == 3)
    {
        if(z_index == 0 || z_index == mExtents[2] - 1)
        {
            return true;
        }
    }
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetExtents(std::vector<unsigned> extents)
{
    if(extents.size()<3)
    {
        EXCEPTION("The extents should be of dimension 3, regardless of element or space dimension");
    }
    mExtents = extents;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetOrigin(c_vector<double, SPACE_DIM> origin)
{
    mOrigin = origin;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetSpacing(double spacing)
{
    mSpacing = spacing;
}

// Explicit instantiation
template class RegularGrid<2> ;
template class RegularGrid<3> ;
