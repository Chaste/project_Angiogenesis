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

#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include "SimpleLinearEllipticSolver.hpp"
#include "VtkMeshWriter.hpp"
#include "ConstBoundaryCondition.hpp"
#include "FiniteElementSolver.hpp"
#include "VesselSurfaceGenerator.hpp"
#include "PlcMesh.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"

#include "KroghCylinderSolver.hpp"

template<unsigned DIM>
KroghCylinderSolver<DIM>::KroghCylinderSolver()
    : AbstractHybridSolver<DIM>(),
      mOuterRadius(100.0),
      mLocations(),
      mpPde(),
      mSolution()
{

}

template<unsigned DIM>
KroghCylinderSolver<DIM>::~KroghCylinderSolver()
{

}

template<unsigned DIM>
std::vector<double> KroghCylinderSolver<DIM>::GetLineSolution()
{
    return mSolution;
}

template<unsigned DIM>
void KroghCylinderSolver<DIM>::SetOuterRadius(double outerRadius)
{
    mOuterRadius = outerRadius;
}

template<unsigned DIM>
void KroghCylinderSolver<DIM>::SetSampleLocations(std::vector<double> locations)
{
    mLocations = locations;
}

template<unsigned DIM>
void KroghCylinderSolver<DIM>::SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
void KroghCylinderSolver<DIM>::Solve(bool writeSolution)
{
    if(!mpPde)
    {
        EXCEPTION("A pde is required to define diffusion and consumption properties.");
    }

    double diffusivity = mpPde->GetDiffusionConstant();

    // If there is a linear in u term assume we are solving the first order problem, otherwise
    // take the zero order version
    double sink_strength;
    double is_first_order = false;
    if(mpPde->GetLinearInUTerm()!=0.0)
    {
        sink_strength = -mpPde->GetLinearInUTerm();
        is_first_order = true;
    }
    else
    {
        sink_strength = -mpPde->GetConstantInUTerm();
    }
    std::string species_label = mpPde->GetVariableName();

    if(!this->mpNetwork)
    {
        EXCEPTION("A vessel network with at least one vessel is required.");
    }
    if(this->mpNetwork->GetNumberOfVessels() < 1)
    {
        EXCEPTION("A vessel network with at least one vessel is required.");
    }
    double capillary_radius = this->mpNetwork->GetVessels()[0]->GetRadius();
    double wall_concentration = this->mpNetwork->GetVessels()[0]->GetSegments()[0]->template GetData<double>(species_label);

    if(!is_first_order)
    {
        double term0 = sink_strength/(2.0 * diffusivity);
        double term1 = 0.0;
        double term2 = 0.0;
        mSolution = std::vector<double>(mLocations.size(), 0.0);
        for(unsigned idx=0; idx<mLocations.size();idx++)
        {
            if(mLocations[idx]>=capillary_radius)
            {
                if(mLocations[idx]<=mOuterRadius)
                {
                    term1 = (term0 / 2.0) * (mLocations[idx]*mLocations[idx] - capillary_radius*capillary_radius);
                    term2 = term0*(mOuterRadius*mOuterRadius)*(std::log(mLocations[idx]/capillary_radius));
                    mSolution[idx] = wall_concentration + term1 - term2;
                }
                else
                {
                    mSolution[idx] = 0.0;
                }
            }
            else
            {
                mSolution[idx] = wall_concentration;
            }
        }
    }
    else
    {
        double term0 = sink_strength /diffusivity;
        double R0 = std::sqrt(term0) * capillary_radius;
        double R1 = std::sqrt(term0) * mOuterRadius;

        mSolution = std::vector<double>(mLocations.size(), 0.0);
        for(unsigned idx=0; idx<mLocations.size();idx++)
        {
            if(mLocations[idx]>=capillary_radius)
            {
                if(mLocations[idx]<=mOuterRadius)
                {
                    double R = std::sqrt(term0) * mLocations[idx];
                    double bessel_k0_R = boost::math::cyl_bessel_k(0, R);
                    double bessel_i0_R = boost::math::cyl_bessel_i(0, R);
                    double bessel_k0_R0 = boost::math::cyl_bessel_k(0, R0);
                    double bessel_i0_R0 = boost::math::cyl_bessel_i(0, R0);
                    double bessel_k1_R1 = boost::math::cyl_bessel_k(1, R1);
                    double bessel_i1_R1 = boost::math::cyl_bessel_i(1, R1);
                    double A = 1.0 * wall_concentration * bessel_k1_R1 / (bessel_i1_R1 * ((bessel_k1_R1 / bessel_i1_R1) * bessel_i0_R0 + bessel_k0_R0));
                    double B = 1.0 * wall_concentration / ((bessel_k1_R1 / bessel_i1_R1) * bessel_i0_R0 + bessel_k0_R0);
                    mSolution[idx] = A * bessel_i0_R + B * bessel_k0_R;
                }
                else
                {
                    mSolution[idx] = 0.0;
                }
            }
            else
            {
                mSolution[idx] = wall_concentration;
            }
        }
    }

    if(writeSolution)
    {
        Write();
    }
}

template<unsigned DIM>
void KroghCylinderSolver<DIM>::Write()
{
    std::ofstream output_file(this->mWorkingDirectory.append("concentration.dat").c_str());
    if (output_file.is_open())
    {
        output_file << "Location Concentration \n";
        for(unsigned idx=0; idx<mLocations.size();idx++)
        {
            output_file << mLocations[idx] << " " << mSolution[idx] << " \n";
        }
        output_file.close();
    }
}

// Explicit instantiation
template class KroghCylinderSolver<2> ;
template class KroghCylinderSolver<3> ;
