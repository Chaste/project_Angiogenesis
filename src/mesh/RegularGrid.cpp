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

RegularGrid::RegularGrid() :
    mSpacing(1.0),
    mExtents(std::vector<unsigned>(3,10)),
    mOrigin(zero_vector<double>(3))
{

}

boost::shared_ptr<RegularGrid> RegularGrid::Create()
{
    MAKE_PTR(RegularGrid, pSelf);
    return pSelf;
}

RegularGrid::~RegularGrid()
{
}

unsigned RegularGrid::Get1dGridIndex(unsigned x_index, unsigned y_index, unsigned z_index)
{
    return x_index + mExtents[0] * y_index + mExtents[0] * mExtents[1] * z_index;
}

std::vector<unsigned> RegularGrid::GetExtents()
{
    return mExtents;
}

c_vector<double, 3> RegularGrid::GetLocation(unsigned x_index, unsigned y_index, unsigned z_index)
{
    c_vector<double, 3> location;
    location[0] = double(x_index) * mSpacing + mOrigin[0];
    location[1] = double(y_index) * mSpacing + mOrigin[1];
    location[2] = double(z_index) * mSpacing + mOrigin[2];
    return location;
}

c_vector<double, 3> RegularGrid::GetLocationOf1dIndex(unsigned grid_index)
{
    unsigned mod_z = grid_index % (mExtents[0] * mExtents[1]);
    unsigned z_index = (grid_index - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;
    return GetLocation(x_index, y_index, z_index);
}

c_vector<double, 3> RegularGrid::GetOrigin()
{
    return mOrigin;
}

unsigned RegularGrid::GetNumberOfPoints()
{
    return mExtents[0] * mExtents[1] * mExtents[2];
}

double RegularGrid::GetSpacing()
{
    return mSpacing;
}


void RegularGrid::SetExtents(std::vector<unsigned> extents)
{
    if(extents.size()<3)
    {
        EXCEPTION("The extents should be of dimension 3");
    }
    mExtents = extents;
}

void RegularGrid::SetOrigin(c_vector<double, 3> origin)
{
    mOrigin = origin;
}

void RegularGrid::SetSpacing(double spacing)
{
    mSpacing = spacing;
}
