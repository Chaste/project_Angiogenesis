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

#include "RandomNumberGenerator.hpp"
#include "GeometryTools.hpp"
#include "OffLatticeRandomNormalSproutingRule.hpp"

template<unsigned DIM>
OffLatticeRandomNormalSproutingRule<DIM>::OffLatticeRandomNormalSproutingRule()
    : AbstractSolutionDependentSproutingRule<DIM>()
{

}

template<unsigned DIM>
OffLatticeRandomNormalSproutingRule<DIM>::~OffLatticeRandomNormalSproutingRule()
{

}

template <unsigned DIM>
boost::shared_ptr<OffLatticeRandomNormalSproutingRule<DIM> > OffLatticeRandomNormalSproutingRule<DIM>::Create()
{
    MAKE_PTR(OffLatticeRandomNormalSproutingRule<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > OffLatticeRandomNormalSproutingRule<DIM>::GetSproutDirection(std::vector<bool> sprout_indices)
{

    std::vector<bool> indices;
    if(sprout_indices.size()>0)
    {
        indices = sprout_indices;
    }
    else
    {
        indices = this->WillSprout();
    }

    std::vector<c_vector<double, DIM> > directions;
    for(unsigned idx = 0; idx < this->mNodes.size(); idx++)
    {
        if(indices[idx])
        {
            c_vector<double, DIM> sprout_direction;
            c_vector<double, DIM> cross_product = VectorProduct(this->mNodes[idx]->GetVesselSegments()[0]->GetUnitTangent(),
                                                                this->mNodes[idx]->GetVesselSegments()[1]->GetUnitTangent());
            double sum = 0.0;
            for(unsigned jdx=0; jdx<DIM; jdx++)
            {
                sum += cross_product[jdx];
            }
            if (sum<=1.e-6)
            {
                // parallel segments, chose any normal to the first tangent
                c_vector<double, DIM> normal;
                c_vector<double, DIM> tangent = this->mNodes[idx]->GetVesselSegments()[0]->GetUnitTangent();

                if(DIM==2 or tangent[2]==0.0)
                {
                    if(tangent[1] == 0.0)
                    {
                        normal[0] = 0.0;
                        normal[1] = 1.0;
                    }
                    else
                    {
                        normal[0] = 1.0;
                        normal[1] = -tangent[0] /tangent[1];
                    }

                }
                else
                {
                    normal[2] = -(tangent[0] + tangent[0])/tangent[2];
                }
                if(RandomNumberGenerator::Instance()->ranf()>=0.5)
                {
                    sprout_direction = normal/norm_2(normal);
                }
                else
                {
                    sprout_direction = -normal/norm_2(normal);
                }
            }
            else
            {
                // otherwise the direction is out of the plane of the segment tangents
                if(RandomNumberGenerator::Instance()->ranf()>=0.5)
                {
                    sprout_direction = cross_product/norm_2(cross_product);
                }
                else
                {
                    sprout_direction = -cross_product/norm_2(cross_product);
                }
            }

            // Rotate by a random angle around the axis
            double angle = RandomNumberGenerator::Instance()->ranf() * 2.0 * M_PI;
            directions.push_back(RotateAboutAxis<DIM>(sprout_direction, this->mNodes[idx]->GetVesselSegments()[0]->GetUnitTangent(), angle));
        }
        else
        {
            directions.push_back(zero_vector<double>(DIM));
        }
    }
    return directions;
}

// Explicit instantiation
template class OffLatticeRandomNormalSproutingRule<2> ;
template class OffLatticeRandomNormalSproutingRule<3> ;
