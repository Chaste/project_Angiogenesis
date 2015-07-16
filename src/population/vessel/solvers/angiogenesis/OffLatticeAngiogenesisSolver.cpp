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

#include <boost/lexical_cast.hpp>
#include "UblasVectorInclude.hpp"
#include "UblasIncludes.hpp"
#include "RandomNumberGenerator.hpp"
#include "VascularNode.hpp"
#include "OffLatticeAngiogenesisSolver.hpp"

template<unsigned DIM>
OffLatticeAngiogenesisSolver<DIM>::OffLatticeAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork, const std::string& rOutputDirectory)
    : AbstractAngiogenesisSolver<DIM>(pNetwork, rOutputDirectory)

{

}

template<unsigned DIM>
OffLatticeAngiogenesisSolver<DIM>::~OffLatticeAngiogenesisSolver()
{

}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeAngiogenesisSolver<DIM>::RotateAboutAxis(c_vector<double, DIM> direction, c_vector<double, DIM> axis, double angle)
{
    double sin_a = std::sin(angle);
    double cos_a = std::cos(angle);
    c_vector<double, DIM> unit_axis = axis / norm_2(axis);

    double dot_product = inner_prod(direction, unit_axis);
    c_vector<double, DIM> new_direction;

    if(DIM==3)
    {
        new_direction[0] = (unit_axis[0] * dot_product * (1.0 - cos_a) + direction[0] * cos_a
                    + (-unit_axis[2] * direction[1] + unit_axis[1] * direction[2]) * sin_a);
        new_direction[1] = (unit_axis[1] * dot_product * (1.0 - cos_a) + direction[1] * cos_a
                    + (unit_axis[2] * direction[0] - unit_axis[0] * direction[2]) * sin_a);
        new_direction[2] = (unit_axis[2] * dot_product * (1.0 - cos_a) + direction[2] * cos_a
                    + (-unit_axis[1] * direction[0] + unit_axis[0] * direction[1]) * sin_a);
    }
    else
    {
        new_direction[0] = unit_axis[0] * dot_product * (1.0 - cos_a) + direction[0] * cos_a;
        new_direction[1] = unit_axis[1] * dot_product * (1.0 - cos_a) + direction[1] * cos_a;
    }
    return new_direction;
}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeAngiogenesisSolver<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection)
{
    c_vector<double, DIM> new_direction;
    c_vector<double, DIM> x_axis = unit_vector<double>(DIM,0);
    c_vector<double, DIM> y_axis = unit_vector<double>(DIM,1);
    c_vector<double, DIM> z_axis;
    if(DIM==3)
    {
        z_axis = unit_vector<double>(DIM,2);
    }

    // Rotate about global axes through random angles
    double angle = M_PI/18.0;
    c_vector<double, DIM> new_directionz = RotateAboutAxis(currentDirection, z_axis, RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, angle));
    c_vector<double, DIM> new_directiony = RotateAboutAxis(new_directionz, y_axis, RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, angle));
    new_direction = RotateAboutAxis(new_directiony, x_axis, RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, angle));

    return new_direction;
}

// Explicit instantiation
template class OffLatticeAngiogenesisSolver<2> ;
template class OffLatticeAngiogenesisSolver<3> ;
