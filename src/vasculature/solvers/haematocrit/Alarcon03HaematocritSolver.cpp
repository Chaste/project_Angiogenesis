///*
//
//Copyright (c) 2005-2015, University of Oxford.
//All rights reserved.
//
//University of Oxford means the Chancellor, Masters and Scholars of the
//University of Oxford, having an administrative office at Wellington
//Square, Oxford OX1 2JD, UK.
//
//This file is part of Chaste.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
//   contributors may be used to endorse or promote products derived from this
//   software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS unsignedERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <stdio.h>
#include "Alarcon03HaematocritSolver.hpp"
#include "LinearSystem.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::Alarcon03HaematocritSolver(double haematocrit) : AbstractHaematocritSolver<DIM>(),
    mTHR(2.5),
    mAlpha(0.5),
    mHaematocrit(haematocrit)
{

}

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::~Alarcon03HaematocritSolver()
{

}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetTHR(double THR)
{
    mTHR = THR;
    assert(mTHR > 1);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetAlpha(double Alpha)
{
    mAlpha = Alpha;
    assert(mAlpha < 1);
    assert(mAlpha > 0);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetHaematocrit(double haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

    // create extra data tables to aid with formation of coefficient matrix for haematocrit calculation
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetVesselEndNodes();
    std::vector< std::vector<unsigned> > VesselsFlowingOutOfNode(nodes.size());
    std::vector< std::vector<unsigned> > VesselsFlowingIntoNode(nodes.size());
    std::vector< std::vector<unsigned> > VesselsAttachedToNodeWithZeroFlow(nodes.size());

    for (unsigned i = 0; i < nodes.size(); i++)
    {
        for (unsigned j = 0; j < nodes[i]->GetNumberOfSegments(); j++)
        {
            if (nodes[i] == nodes[i]->GetVesselSegment(j)->GetVessel()->GetStartNode())
            {
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() < 0)
                {
                    VesselsFlowingIntoNode[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() > 0)
                {
                    VesselsFlowingOutOfNode[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() == 0)
                {
                    VesselsAttachedToNodeWithZeroFlow[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
            }
            if (nodes[i] == nodes[i]->GetVesselSegment(j)->GetVessel()->GetEndNode() && nodes[i]->GetVesselSegment(j)->GetVessel()->GetEndNode() !=
                    nodes[i]->GetVesselSegment(j)->GetVessel()->GetStartNode())
            {
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() > 0)
                {
                    VesselsFlowingIntoNode[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() < 0)
                {
                    VesselsFlowingOutOfNode[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
                if (nodes[i]->GetVesselSegment(j)->GetFlowVelocity() == 0)
                {
                    VesselsAttachedToNodeWithZeroFlow[i].push_back(vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(j)->GetVessel()));
                }
            }
        }

        assert(nodes[i]->GetNumberOfSegments() == VesselsAttachedToNodeWithZeroFlow[i].size() + VesselsFlowingOutOfNode[i].size() + VesselsFlowingIntoNode[i].size());
    }

    // Set up the system
    PetscInt lhsVectorSize = vascularNetwork->GetNumberOfVessels();
    LinearSystem linearSystem(lhsVectorSize, 3);
    linearSystem.SetPcType("lu");
    PetscOptionsSetValue("-pc_factor_mat_solver_package", "umfpack");
    PetscOptionsSetValue("-pc_factor_zeropivot", 0);
    linearSystem.SetKspType("gmres");

    unsigned EquationNumber = 0;

    // equations which say that arterial input vessels have an arterial haematocrit level

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVesselNodes(); i++)
    {
        if (nodes[i]->IsInputNode())
        {
            linearSystem.AssembleIntermediateLinearSystem();
            linearSystem.AddToMatrixElement(EquationNumber, vascularNetwork->GetVesselIndex(nodes[i]->GetVesselSegment(0)->GetVessel()), 1);
            linearSystem.AssembleIntermediateLinearSystem();
            linearSystem.SetRhsVectorElement(EquationNumber, mHaematocrit);
            EquationNumber++;
        }
    }

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVesselNodes(); i++)
    {

        if (VesselsFlowingIntoNode[i].size() + VesselsFlowingOutOfNode[i].size() + VesselsAttachedToNodeWithZeroFlow[i].size() > 3)
        {
            unsigned totalNumberOfVesselsAttachedToNode = VesselsFlowingIntoNode[i].size() + VesselsFlowingOutOfNode[i].size() + VesselsAttachedToNodeWithZeroFlow[i].size();

            cout << "Error: Alarcon03HaematocritSolver is not equipped to handle vessel networks where more than 3 vessels unsignedersect at a node pounsigned in the network. Perhaps replace this calculation with a constant haematocrit calculation.\n";
            assert(totalNumberOfVesselsAttachedToNode <= 3);
        }

        if (VesselsFlowingIntoNode[i].size() + VesselsFlowingOutOfNode[i].size() == 2)
        {
            // must be one flow going in to node and one flowing out for conservation
            assert(VesselsFlowingIntoNode[i].size() == 1 && VesselsFlowingOutOfNode[i].size() == 1);

            linearSystem.AssembleIntermediateLinearSystem();
            linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][0], 1);
            linearSystem.AssembleIntermediateLinearSystem();
            linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], -1);
            EquationNumber++;
        }

        if (VesselsFlowingIntoNode[i].size() + VesselsFlowingOutOfNode[i].size() == 3)
        {

            if (VesselsFlowingIntoNode[i].size() == 2)
            {
                // conservation equation for node
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][0], 1);
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][1], 1);
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], -1);
                EquationNumber++;

            }
            if (VesselsFlowingIntoNode[i].size() == 1)
            {
                // conservation equation for node
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][0], 1);
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], -1);
                linearSystem.AssembleIntermediateLinearSystem();
                linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][1], -1);
                EquationNumber++;


                if (fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity()) >= fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity()))
                {
                    if (fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity()) < mTHR*fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity()))
                    {
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], 1);
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][1], -mAlpha*(fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity())/fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity())));
                    }
                    else
                    {
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], 1);
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][0], -1);
                    }

                }

                if (fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity()) > fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity()))
                {
                    if (fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity()) < mTHR*fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity()))
                    {
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][1], 1);
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][0], -mAlpha*(fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowVelocity())/fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowVelocity())));
                    }
                    else
                    {
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingOutOfNode[i][1], 1);
                        linearSystem.AssembleIntermediateLinearSystem();
                        linearSystem.AddToMatrixElement(EquationNumber, VesselsFlowingIntoNode[i][0], -1);
                    }

                }


                EquationNumber++;

            }

        }

    }

    // equations which say that vessels with zero flow have zero haematocrit in

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVessels(); i++)
    {
        if (vascularNetwork->GetVessel(i)->GetFlowVelocity() == 0)
        {
            linearSystem.AssembleIntermediateLinearSystem();
            linearSystem.SetMatrixElement(EquationNumber, i, 1.0);
            //            linearSystem.AssembleIntermediateLinearSystem();
            //            linearSystem.SetRhsVectorElement(EquationNumber, 0.0);
            EquationNumber++;
        }
    }

    linearSystem.AssembleFinalLinearSystem();

    Vec solution = PetscTools::CreateVec(vascularNetwork->GetNumberOfVessels());
    Vec initialGuess = PetscTools::CreateVec(vascularNetwork->GetNumberOfVessels());

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVessels(); i++)
    {
        if (vascularNetwork->GetVessel(i)->GetFlowRate() == 0)
        {
            PetscVecTools::SetElement(initialGuess, i, 0);
        }
        else
        {
            PetscVecTools::SetElement(initialGuess, i, vascularNetwork->GetVessel(i)->GetHaematocrit());
        }
    }

    try
    {
        solution = linearSystem.Solve(initialGuess);
    }
    catch (Exception &e)
    {
        std::cout << e.GetMessage() << std::endl;
    }

    // deal with minor rounding errors in calculation

    ReplicatableVector a(solution);

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVessels(); i++)
    {
        if (a[i] < pow(10.0,-15))
        {
            a[i] = 0;
        }

    }

    // assign haematocrit levels to vessels

    for (unsigned i = 0; i < vascularNetwork->GetNumberOfVessels(); i++)
    {
        for (unsigned jdx = 0; jdx < vascularNetwork->GetVessel(i)->GetNumberOfSegments(); jdx++)
        {
            vascularNetwork->GetVessel(i)->GetSegment(jdx)->SetHaematocrit(double(a[i]));
        }
    }

    /*
     * clean up
     */
    PetscTools::Destroy(solution);
    PetscTools::Destroy(initialGuess);

}


// Explicit instantiation

template class Alarcon03HaematocritSolver<2>;
template class Alarcon03HaematocritSolver<3>;
