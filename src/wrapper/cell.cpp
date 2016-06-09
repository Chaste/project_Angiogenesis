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

#ifdef CHASTE_ANGIOGENESIS_PYTHON
#include <vector>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "SmartPointers.hpp"
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"
#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"
#include "CellsGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellId.hpp"
#include "UblasIncludes.hpp"
#include "SimpleCell.hpp"
#include "SimpleCellPopulation.hpp"
#include "CancerCellMutationState.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractMesh.hpp"
#include "AbstractCellPopulation.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
//#include "TestCases.cpp"

#include "converters.hpp"

using namespace boost::python;

// Cell Overloads
const c_vector<double, 3>& (SimpleCell<3>::*SimpleCell_rGetLocation_3)() const = &SimpleCell<3>::rGetLocation;

boost::shared_ptr<Cell> GenerateCell(const std::string& mutation_state, const std::string& cell_cycle_model)
{
    MAKE_PTR(CancerCellMutationState, p_state); // Default state

    Owen2011OxygenBasedCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;
    p_model->SetDimension(2);

    CellPtr p_cell(new Cell(p_state, p_model));
    p_cell->SetApoptosisTime(30.0);
    p_cell->GetCellData()->SetItem("oxygen", 30.0);

    MAKE_PTR(StemCellProliferativeType, p_stem_type); // Default Type
    p_cell->SetCellProliferativeType(p_stem_type);
    return p_cell;
}

boost::shared_ptr<CaBasedCellPopulation<3> > GenerateCaBasedCellPopulation(boost::shared_ptr<PottsMesh<3u> > pMesh, std::vector<CellPtr> cells, std::vector<unsigned> location_indices)
{
    boost::shared_ptr<CaBasedCellPopulation<3> > p_cell_population(new CaBasedCellPopulation<3>(*pMesh, cells, location_indices));
    return p_cell_population;
}


// Make the module
BOOST_PYTHON_MODULE(_cell)
{
    class_<SimpleCell<3>, boost::shared_ptr<SimpleCell<3> > >("SimpleCell", init<double, double, double>())
       .def(init<c_vector<double, 3> >())
       .def("GetLocation", SimpleCell_rGetLocation_3, return_value_policy<copy_const_reference>())
   ;

    class_<SimpleCellPopulation<3>, boost::shared_ptr<SimpleCellPopulation<3> > >("SimpleCellPopulation")
        .def("AddCells", &SimpleCellPopulation<3>::AddCells)
        .def("AddCell", &SimpleCellPopulation<3>::AddCell)
        .def("GetCells", &SimpleCellPopulation<3>::GetCells)
        .def("Write", &SimpleCellPopulation<3>::Write)
    ;

    class_<Cell, boost::shared_ptr<Cell>, boost::noncopyable>("Cell", init<boost::shared_ptr<AbstractCellProperty>,
            AbstractCellCycleModel*, AbstractSrnModel*, bool, CellPropertyCollection>())
            .def("GetCellId", &Cell::GetCellId)
            .def("SetCellProliferativeType", &Cell::SetCellProliferativeType)
            .def("SetApoptosisTime", &Cell::SetApoptosisTime)
            .def("GetCellData", &Cell::GetCellData)
    ;

    class_<CaBasedCellPopulation<3>, boost::shared_ptr<CaBasedCellPopulation<3> >, boost::noncopyable>("CaBasedCellPopulation",
            init<PottsMesh<3>&, std::vector<CellPtr>&, const std::vector<unsigned>, unsigned, bool, bool>())
            .def("GetNumNodes", &CaBasedCellPopulation<3>::GetNumNodes)
            .def("GetNumRealCells", &CaBasedCellPopulation<3>::GetNumRealCells)
    ;

    def("CellGenerator", &GenerateCell);

    def("CaBasedCellPopulationGenerator", &GenerateCaBasedCellPopulation);

    class_<CancerCellMutationState, boost::shared_ptr<CancerCellMutationState>, boost::noncopyable>("CancerCellMutationState", init<>())
    ;

    class_<std::vector<boost::shared_ptr<Cell> > > ("VecCellPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<Cell > > >())
    ;

    class_<std::vector<boost::shared_ptr<SimpleCell<3> > > > ("VecSimpleCellPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<SimpleCell<3> > > >())
   ;

    class_<std::vector<boost::shared_ptr<SimpleCellPopulation<3> > > > ("VecSimpleCellPopulationPtrs")
         .def(vector_ptr_indexing_suite<std::vector<boost::shared_ptr<SimpleCellPopulation<3> > > >())
    ;

    // Register Angiogenesis Converters
    PythonIterableToStl()
        .from_python<std::vector<boost::shared_ptr<Cell > > >()
        .from_python<std::vector<boost::shared_ptr<SimpleCell<3> > > >()
      ;
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
