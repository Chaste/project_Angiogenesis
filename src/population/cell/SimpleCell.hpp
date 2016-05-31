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

#ifndef SIMPLECELL_HPP_
#define SIMPLECELL_HPP_

#include "ChastePoint.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellProperty.hpp"

/* A minimal cell class. Can be used to create Chaste cells of different types.
 */

template<unsigned DIM>
class SimpleCell : public ChastePoint<DIM>
{
    unsigned mIndex;

    std::string mCellCycleModel;

    std::string mCellMutationState;

    std::string mCellProliferativeType;

public:

    /* Constructor
     */
    SimpleCell(double v1 = 0, double v2 = 0, double v3 = 0);

    /* Constructor
     */
    SimpleCell(c_vector<double, DIM> location);

    /* Factory constructor method
     * @return a shared pointer to a new cell
     */
    static boost::shared_ptr<SimpleCell<DIM> > Create(double v1 = 0, double v2 = 0, double v3 = 0);

    /* Factory constructor method
     * @return a shared pointer to a new cell
     */
    static boost::shared_ptr<SimpleCell<DIM> > Create(c_vector<double, DIM> location);

    /* Desctructor
     */
    ~SimpleCell();

    /* Return the index
     */
    unsigned GetIndex();

    /* Set the index
     */
    void SetIndex(unsigned index);

    const std::string& GetCellCycleModel();

    const std::string& GetCellMutationState();

    const std::string& GetCellProliferativeType();

    void SetCellCycleModel(std::string& modelName);

    void SetCellMutationState(std::string& mutationState);

    void SetCellProliferativeType(std::string& proliferativeType);
};

#endif /* SIMPLECELL_HPP_*/
