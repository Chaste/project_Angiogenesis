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

#include "../growth_direction_modifiers/OffLatticePrwGrowthDirectionModifier.hpp"

#include "GeometryTools.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
OffLatticePrwGrowthDirectionModifier<DIM>::OffLatticePrwGrowthDirectionModifier()
    : AbstractGrowthDirectionModifier<DIM>(),
      mGlobalX(unit_vector<double>(DIM,0)),
      mGlobalY(unit_vector<double>(DIM,0)),
      mGlobalZ(zero_vector<double>(DIM)),
      mMeanAngles(std::vector<double>(DIM, 0.0)),
      mSdvAngles(std::vector<double>(DIM, M_PI/18.0))
{
    if(DIM==3)
    {
        mGlobalZ = unit_vector<double>(DIM,2);
    }
}

template <unsigned DIM>
boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<DIM> > OffLatticePrwGrowthDirectionModifier<DIM>::Create()
{
    MAKE_PTR(OffLatticePrwGrowthDirectionModifier<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
OffLatticePrwGrowthDirectionModifier<DIM>::~OffLatticePrwGrowthDirectionModifier()
{

}

template<unsigned DIM>
c_vector<double, DIM> OffLatticePrwGrowthDirectionModifier<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                                                                    boost::shared_ptr<VascularNode<DIM> > pNode)
{
    double angle_x = RandomNumberGenerator::Instance()->NormalRandomDeviate(mMeanAngles[0], mSdvAngles[0]);
    double angle_y = RandomNumberGenerator::Instance()->NormalRandomDeviate(mMeanAngles[1], mSdvAngles[1]);
    double angle_z = 0.0;
    if(DIM==3)
    {
        angle_z = RandomNumberGenerator::Instance()->NormalRandomDeviate(mMeanAngles[2], mSdvAngles[2]);
    }

    c_vector<double, DIM> new_direction_z = RotateAboutAxis<DIM>(currentDirection, mGlobalZ, angle_z);
    c_vector<double, DIM> new_direction_y = RotateAboutAxis<DIM>(new_direction_z, mGlobalY, angle_y);
    return RotateAboutAxis<DIM>(new_direction_y, mGlobalX, angle_x);
}

// Explicit instantiation
template class OffLatticePrwGrowthDirectionModifier<2> ;
template class OffLatticePrwGrowthDirectionModifier<3> ;
