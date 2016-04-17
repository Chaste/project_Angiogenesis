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

#ifndef _OWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define _OWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_


#include <cmath>

/**
 * Represents the Owen (2011) system of ODEs
 *
 * The variables are
 * 0. phi = Cell-cycle phase
 * 1. p53 = p53 concentration
 * 2. VEGF = VEGF concentration
 * 3. O2 = Oxygen concentration
 */

class Owen2011OxygenBasedCellCycleOdeSystem
{
private:


public:

    /**
     * Constructor.
     *
     * @param oxygenConcentration is a non-dimensional oxygen concentration value between 0 and 1
     * @param isLabelled whether the cell associated with this cell cycle ODE system is labelled (this affects the ODE system)
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    Owen2011OxygenBasedCellCycleOdeSystem();

    /**
     * Destructor.
     */
    ~Owen2011OxygenBasedCellCycleOdeSystem();

    
};

// Declare identifier for the serializer
//#include "SerializationExportWrapper.hpp"
//CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleOdeSystem)


#endif /*_OWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
