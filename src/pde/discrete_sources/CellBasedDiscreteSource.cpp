/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "CellBasedDiscreteSource.hpp"
#include "AbstractCellPopulation.hpp"
#include "VesselNetwork.hpp"
#include "GeometryTools.hpp"

template<unsigned DIM>
CellBasedDiscreteSource<DIM>::CellBasedDiscreteSource()
    :   DiscreteSource<DIM>(),
        mCellConstantInUValue(0.0*unit::mole_per_second),
        mCellLinearInUValue(0.0*unit::per_second)
{

}

template<unsigned DIM>
CellBasedDiscreteSource<DIM>::~CellBasedDiscreteSource()
{

}

template<unsigned DIM>
boost::shared_ptr<CellBasedDiscreteSource<DIM> > CellBasedDiscreteSource<DIM>::Create()
{
    MAKE_PTR(CellBasedDiscreteSource<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
std::vector<units::quantity<unit::concentration_flow_rate> > CellBasedDiscreteSource<DIM>::GetConstantInUMeshValues()
{
    if(!this->mpMesh)
    {
        EXCEPTION("A mesh is required for this type of source");
    }
    return std::vector<units::quantity<unit::concentration_flow_rate> >();
}

template<unsigned DIM>
std::vector<units::quantity<unit::rate> > CellBasedDiscreteSource<DIM>::GetLinearInUMeshValues()
{
    if(!this->mpMesh)
    {
        EXCEPTION("A mesh is required for this type of source");
    }
    return std::vector<units::quantity<unit::rate> >();
}

template<unsigned DIM>
std::vector<units::quantity<unit::concentration_flow_rate> > CellBasedDiscreteSource<DIM>::GetConstantInURegularGridValues()
{
    if(!this->mpRegularGrid)
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    std::vector<units::quantity<unit::concentration_flow_rate> > values(this->mpRegularGrid->GetNumberOfPoints(), 0.0*unit::mole_per_metre_cubed_per_second);
    units::quantity<unit::length> grid_spacing = this->mpRegularGrid->GetSpacing();
    units::quantity<unit::volume> grid_volume = units::pow<3>(grid_spacing);

    std::vector<std::vector<CellPtr> > point_cell_map = this->mpRegularGrid->GetPointCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        values[idx] += mCellConstantInUValue * double(point_cell_map[idx].size())/grid_volume;
    }
    return values;

}

template<unsigned DIM>
std::vector<units::quantity<unit::rate> > CellBasedDiscreteSource<DIM>::GetLinearInURegularGridValues()
{
    if(!this->mpRegularGrid)
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    std::vector<units::quantity<unit::rate> > values(this->mpRegularGrid->GetNumberOfPoints(), 0.0*unit::per_second);
    std::vector<std::vector<CellPtr> > point_cell_map = this->mpRegularGrid->GetPointCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        values[idx] += mCellLinearInUValue * double(point_cell_map[idx].size());
    }
    return values;
}

template<unsigned DIM>
void CellBasedDiscreteSource<DIM>::SetConstantInUConsumptionRatePerCell(units::quantity<unit::molar_flow_rate> value)
{
    mCellConstantInUValue = value;
}

template<unsigned DIM>
void CellBasedDiscreteSource<DIM>::SetLinearInUConsumptionRatePerCell(units::quantity<unit::rate> value)
{
    mCellLinearInUValue = value;
}

// Explicit instantiation
template class CellBasedDiscreteSource<2>;
template class CellBasedDiscreteSource<3>;
