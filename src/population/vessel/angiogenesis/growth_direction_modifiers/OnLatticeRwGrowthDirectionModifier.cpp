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

#include "../growth_direction_modifiers/OnLatticeRwGrowthDirectionModifier.hpp"

#include "GeometryTools.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
OnLatticeRwGrowthDirectionModifier<DIM>::OnLatticeRwGrowthDirectionModifier()
    : AbstractGrowthDirectionModifier<DIM>(),
      mGlobalX(unit_vector<double>(DIM,0)),
      mGlobalY(unit_vector<double>(DIM,0)),
      mGlobalZ(zero_vector<double>(DIM))

{
    if(DIM==3)
    {
        mGlobalZ = unit_vector<double>(DIM,2);
    }
}

template<unsigned DIM>
OnLatticeRwGrowthDirectionModifier<DIM>::~OnLatticeRwGrowthDirectionModifier()
{

}

template<unsigned DIM>
c_vector<double, DIM> OnLatticeRwGrowthDirectionModifier<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                                                                          boost::shared_ptr<VascularNode<DIM> > pNode)
{
    // Assign an index to each of the possible direction
    unsigned direction_index = 0;
    if(-currentDirection[0] > 0 && -currentDirection[0] > 1.0 -1.e-6)
    {
        direction_index = 1;
    }
    else if(-currentDirection[0] < 0 &&  -currentDirection[0] < -1.0 +1.e-6)
    {
        direction_index = 2;
    }
    else if(-currentDirection[1] > 0 && -currentDirection[1] > 1.0 -1.e-6)
    {
        direction_index = 3;
    }
    else if(-currentDirection[1] < 0 &&  -currentDirection[1] < -1.0 +1.e-6)
    {
        direction_index = 4;
    }
    if(DIM ==3)
    {
        if(-currentDirection[2] > 0 && -currentDirection[2] > 1.0 -1.e-6)
        {
            direction_index = 5;
        }
        if(-currentDirection[2] < 0 &&  -currentDirection[2] < -1.0 +1.e-6)
        {
            direction_index = 6;
        }
    }
    if(direction_index == 0)
    {
        EXCEPTION("Initial network should be aligned with the grid.");
    }

    // Get a new direction that does not go back on itself
    bool found = false;
    unsigned new_direction_index = 0;
    while(!found)
    {
        unsigned random_num = (RandomNumberGenerator::Instance()->randMod(2*DIM) + 1);
        if(random_num != direction_index)
        {
            new_direction_index = random_num;
            found = true;
            break;
        }
    }

    // Get the direction from the index
    c_vector<double, DIM> new_direction;
    if(new_direction_index == 1)
    {
        new_direction = unit_vector<double>(DIM,0);
    }
    else if(new_direction_index == 2)
    {
        new_direction = -unit_vector<double>(DIM,0);
    }
    else if(new_direction_index == 3)
    {
        new_direction = unit_vector<double>(DIM,1);
    }
    else if(new_direction_index == 4)
    {
        new_direction = -unit_vector<double>(DIM,1);
    }

    if(DIM==3)
    {
        if(new_direction_index == 5)
        {
            new_direction = unit_vector<double>(DIM,2);
        }
        else if(new_direction_index == 6)
        {
            new_direction = -unit_vector<double>(DIM,2);
        }
    }
    return new_direction;
}

// Explicit instantiation
template class OnLatticeRwGrowthDirectionModifier<2> ;
template class OnLatticeRwGrowthDirectionModifier<3> ;
