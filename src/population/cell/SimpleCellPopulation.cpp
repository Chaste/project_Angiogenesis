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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyData.h"

#include "SimpleCellPopulation.hpp"

SimpleCellPopulation::SimpleCellPopulation() :
        mCells()
{

}

boost::shared_ptr<SimpleCellPopulation> SimpleCellPopulation::Create()
{
    MAKE_PTR(SimpleCellPopulation, pSelf);
    return pSelf;
}

SimpleCellPopulation::~SimpleCellPopulation()
{

}

class ExistsInVector
{

    std::vector<boost::shared_ptr<SimpleCell> > m_vec;

public:

    ExistsInVector(std::vector<boost::shared_ptr<SimpleCell> > vec)
        : m_vec(vec)
    {

    }
    bool operator() (boost::shared_ptr<SimpleCell> i)
    {
        return (std::find(m_vec.begin(), m_vec.end(), i) != m_vec.end());
    }
};

void SimpleCellPopulation::BooleanWithVesselNetwork(boost::shared_ptr<CaVascularNetwork<3> > pNetwork)
{
    std::vector<boost::shared_ptr<SimpleCell> > remove_cells;
    double tolerance = 1.e-6;
    for(unsigned idx=0; idx<mCells.size();idx++)
    {
        std::pair<boost::shared_ptr<CaVesselSegment<3> >, double> seg_pair = pNetwork->GetNearestSegment(mCells[idx]->rGetLocation());
        if(seg_pair.second < tolerance)
        {
            remove_cells.push_back(mCells[idx]);
        }
    }

    mCells.erase( std::remove_if(mCells.begin(), mCells.end(), ExistsInVector(remove_cells)), mCells.end());
}

void SimpleCellPopulation::GenerateCellsOnGrid(unsigned xDim, unsigned yDim, unsigned zDim, double spacing, c_vector<double, 3> origin)
{
    std::vector<boost::shared_ptr<SimpleCell> > cells;
    for (unsigned idx = 0; idx < zDim;idx++)
    {
        for (unsigned jdx = 0; jdx < yDim ;jdx++)
        {
            for (unsigned kdx = 0; kdx < xDim;kdx++)
            {
                double x_coord = double(kdx) * spacing + origin[0];
                double y_coord = double(jdx) * spacing + origin[1];
                double z_coord = double(idx) * spacing + origin[2];
                cells.push_back(SimpleCell::Create(x_coord, y_coord, z_coord));
            }
        }
    }
    AddCells(cells);
}

void SimpleCellPopulation::GenerateCellsOnGrid(boost::shared_ptr<Part> pPart, double spacing)
{
    std::vector<boost::shared_ptr<SimpleCell> > cells;

    // Get the bounding box of the part
    c_vector<double,6> bbox = pPart->GetBoundingBox();
    unsigned num_x = double(bbox[1] - bbox[0]) / spacing + 1;
    unsigned num_y = double(bbox[3] - bbox[2]) / spacing + 1;
    unsigned num_z = double(bbox[5] - bbox[4]) / spacing + 1;

    bool first_loop = true;
    for (unsigned idx = 0; idx < num_z;idx++)
    {
        for (unsigned jdx = 0; jdx < num_y ;jdx++)
        {
            for (unsigned kdx = 0; kdx < num_x;kdx++)
            {
                c_vector<double, 3> location;
                location[0] = bbox[0] + double(kdx) * spacing;
                location[1] = bbox[2] + double(jdx) * spacing;
                location[2] = bbox[4] + double(idx) * spacing;

                if(pPart->IsPointInPart(location, first_loop))
                {
                    cells.push_back(SimpleCell::Create(location));
                }
                first_loop = false;
            }
        }
    }
    AddCells(cells);
}

std::vector<boost::shared_ptr<SimpleCell> > SimpleCellPopulation::GetCells()
{
    return mCells;
}

void SimpleCellPopulation::AddCell(boost::shared_ptr<SimpleCell> pCell)
{
    mCells.push_back(pCell);
}

void SimpleCellPopulation::AddCells(std::vector<boost::shared_ptr<SimpleCell> > cells)
{
    mCells.insert(mCells.end(), cells.begin(), cells.end());
}

vtkSmartPointer<vtkPoints> SimpleCellPopulation::GetVtk()
{
    vtkSmartPointer<vtkPoints> p_vertices = vtkSmartPointer<vtkPoints>::New();

    p_vertices->SetNumberOfPoints(mCells.size());
    for (vtkIdType idx = 0; idx < vtkIdType(mCells.size()); idx++)
    {
        c_vector<double, 3> location = mCells[idx]->rGetLocation();
        p_vertices->SetPoint(idx, location[0], location[1], location[2]);
    }

    return p_vertices;
}

void SimpleCellPopulation::Write(const std::string& rFileName)
{
    vtkSmartPointer<vtkPolyData> p_part_data = vtkSmartPointer<vtkPolyData>::New();
    p_part_data->SetPoints(GetVtk());
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(rFileName.c_str());
    writer->SetInput(p_part_data);
    writer->Write();
}
