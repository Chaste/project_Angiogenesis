///*
//
//Copyright (c) 2005-2014, University of Oxford.
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
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//*/
//
//#ifndef ABSTRACTVASCULARSIMULATION_HPP_
//#define ABSTRACTVASCULARSIMULATION_HPP_
//
//#include "AbstractCellPopulation.hpp"
//#include "SimpleCellPopulation.hpp"
//#include "SmartPointers.hpp"
//
///*
// * Class for simulations involving vessels. More generic than the Cell Based
// * simulation classes, as Cell Based Chaste functionality is not always used
// * for vessel based simulations.
// */
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//class AbstractVascularSimulation
//{
//
//    boost::shared_ptr<AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> > mpCellBasedPopulation;
//    boost::shared_ptr<SimpleCellPopulation<SPACE_DIM> > mpSimpleCellPopulation;
//
//public:
//
//    AbstractVascularSimulation();
//
//    virtual ~AbstractVascularSimulation();
//
//    void SetCellBasedPopulation(boost::shared_ptr<AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> > pCellPopulation);
//
//    void SetSimpleCellPopulation(boost::shared_ptr<SimpleCellPopulation<SPACE_DIM> > pCellPopulation);
//
//    void SetCellBasedPdeHandler()
//
//
//
//};
//
//#endif /*ABSTRACTVASCULARSIMULATION_HPP_ */
