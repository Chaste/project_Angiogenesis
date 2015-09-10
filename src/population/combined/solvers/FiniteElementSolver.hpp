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

#ifndef FINITEELEMENTSOLVER_HPP_
#define FINITEELEMENTSOLVER_HPP_

#include "SmartPointers.hpp"
#include "AbstractHybridSolver.hpp"
#include "PlcMesh.hpp"
#include "Part.hpp"
#include "HybridLinearEllipticPde.hpp"

template<unsigned DIM>
class FiniteElementSolver : public AbstractHybridSolver<DIM>
{
    using AbstractHybridSolver<DIM>::Solve;
    boost::shared_ptr<Part<DIM> > mpDomain;
    double mGridSize;
    boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > mpPde;

public:

    FiniteElementSolver();

    ~FiniteElementSolver();

    void SetDomain(boost::shared_ptr<Part<DIM> > pDomain);

    void SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde);

    void SetMaxElementArea(double maxElementArea);

    void Solve(bool writeSolution = false);

private:

    void Write(std::vector<double> output, boost::shared_ptr<PlcMesh<DIM, DIM> > p_mesh);
};

#endif /* FINITEELEMENTSOLVER_HPP_ */
