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

#ifdef CHASTE_VTK
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#endif // CHASTE_VTK
#include "SmartPointers.hpp"
#include "Exception.hpp"

template <unsigned DIM>
VtkVesselNetworkWriter<DIM>::VtkVesselNetworkWriter() :
    mpVesselNetwork(),
    mpVtkVesselNetwork(vtkSmartPointer<vtkPolyData>::New()),
    mIsVtkNetworkUpToDate(false),
    mFilename()
{

}

template <unsigned DIM>
VtkVesselNetworkWriter<DIM>::~VtkVesselNetworkWriter()
{

}

template <unsigned DIM>
boost::shared_ptr<VtkVesselNetworkWriter<DIM> > VtkVesselNetworkWriter<DIM>::Create()
{
    MAKE_PTR(VtkVesselNetworkWriter<DIM>, pSelf);
    return pSelf;
}

template <unsigned DIM>
void VtkVesselNetworkWriter<DIM>::SetVesselNetwork(boost::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    mpVesselNetwork = pNetwork;
    mIsVtkNetworkUpToDate = false;
}

template <unsigned DIM>
vtkSmartPointer<vtkPolyData> VtkVesselNetworkWriter<DIM>::GetOutput()
{
    if(!mpVesselNetwork)
    {
        EXCEPTION("A vessel network is required for the vtk writer.");
    }

    mpVtkVesselNetwork = vtkSmartPointer<vtkPolyData>::New();

    if(mpVesselNetwork->GetNumberOfVessels()>0)
    {
        // Set up the vessel and node data arrays.
        std::vector<vtkSmartPointer<vtkDoubleArray> > pVesselInfoVector;
        std::map<std::string, double>::iterator vessel_map_iterator;
        std::map<std::string, double> vessel_data_map = mpVesselNetwork->GetVessels()[0]->GetOutputData();
        for(vessel_map_iterator = vessel_data_map.begin(); vessel_map_iterator != vessel_data_map.end(); vessel_map_iterator++)
        {
            vtkSmartPointer<vtkDoubleArray> pVesselInfo = vtkSmartPointer<vtkDoubleArray>::New();
            pVesselInfo->SetNumberOfComponents(1);
            pVesselInfo->SetNumberOfTuples(mpVesselNetwork->GetNumberOfVessels());
            pVesselInfo->SetName((*vessel_map_iterator).first.c_str());
            pVesselInfoVector.push_back(pVesselInfo);
        }

        // Set up node info arrays
        std::vector<vtkSmartPointer<vtkDoubleArray> > pNodeInfoVector;
        unsigned numberOfNodes = mpVesselNetwork->GetNumberOfNodes();

        std::map<std::string, double>::iterator vtk_node_map_iterator;
        std::map<std::string, double> vtk_node_map = mpVesselNetwork->GetVessels()[0]->GetStartNode()->GetOutputData();
        for(vtk_node_map_iterator = vtk_node_map.begin(); vtk_node_map_iterator != vtk_node_map.end(); vtk_node_map_iterator++)
        {
            vtkSmartPointer<vtkDoubleArray> pNodeInfo = vtkSmartPointer<vtkDoubleArray>::New();
            pNodeInfo->SetNumberOfComponents(1);
            pNodeInfo->SetNumberOfTuples(numberOfNodes);
            pNodeInfo->SetName((*vtk_node_map_iterator).first.c_str());
            pNodeInfoVector.push_back(pNodeInfo);
        }

        // Create the geometric data
        vtkSmartPointer<vtkPoints> pPoints= vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> pLines = vtkSmartPointer<vtkCellArray>::New();

        unsigned vessel_index=0;
        std::vector<boost::shared_ptr<VesselNode<DIM> > > nodes = this->GetNodes();
        for(unsigned idx=0; idx<nodes.size(); idx++)
        {
            nodes[idx]->SetId(idx);
            if(DIM == 2)
            {
                pPoints->InsertNextPoint(nodes[idx]->rGetLocation()[0], nodes[idx]->rGetLocation()[1], 0.0);
            }
            else
            {
                pPoints->InsertNextPoint(nodes[idx]->rGetLocation()[0], nodes[idx]->rGetLocation()[1], nodes[idx]->rGetLocation()[2]);
            }

            std::map<std::string, double> vtk_node_data = nodes[idx]->GetOutputData();
    //        std::map<std::string, boost::any> generic_node_data = nodes[idx]->rGetDataContainer().GetMap();

            // Add the node data
            for(unsigned jdx=0; jdx < pNodeInfoVector.size(); jdx++)
            {
                // Get the key
                std::string key = pNodeInfoVector[jdx]->GetName();

                // If it is in the vtk data use it
                if(vtk_node_data.count(key) == 1)
                {
                    pNodeInfoVector[jdx]->SetValue(idx, vtk_node_data[key]);
                }
                // Otherwise check the generic data
    //            else if(generic_node_data.count(key) == 1)
    //            {
    //                if(generic_node_data[key].type() == typeid(double))
    //                {
    //                    double cast_value = boost::any_cast<double>(generic_node_data[key]);
    //                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
    //                }
    //                else if(generic_node_data[key].type() == typeid(unsigned))
    //                {
    //                    double cast_value = double(boost::any_cast<unsigned>(generic_node_data[key]));
    //                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
    //                }
    //                else if(generic_node_data[key].type() == typeid(bool))
    //                {
    //                    double cast_value = double(boost::any_cast<bool>(generic_node_data[key]));
    //                    pNodeInfoVector[jdx]->SetValue(idx, cast_value);
    //                }
    //            }
            }
        }

        // Do the vessels
        typename std::vector<boost::shared_ptr<Vessel<DIM> > >::iterator it;
        for(it = mVessels.begin(); it < mVessels.end(); it++)
        {
            vtkSmartPointer<vtkLine> pLine = vtkSmartPointer<vtkLine>::New();
            std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = (*it)->GetSegments();

            for(unsigned i = 0; i < segments.size(); i++)
            {
                pLine->GetPointIds()->InsertId(i, segments[i]->GetNode(0)->GetId());

                // Do an extra insert for the last node in the segment
                if (i == segments.size() - 1)
                {
                    pLine->GetPointIds()->InsertId(i + 1, segments[i]->GetNode(1)->GetId());
                }
            }
            pLines->InsertNextCell(pLine);

            // Add the vessel data
            std::map<std::string, double> vtk_vessel_data = (*it)->GetOutputData();
    //        std::map<std::string, double> generic_vessel_data = (*it)->rGetDataContainer().GetMap();

            for(unsigned idx=0; idx < pVesselInfoVector.size(); idx++)
            {
                // Get the key
                std::string key = pVesselInfoVector[idx]->GetName();
                // If it is in the vtk data use it
                if(vtk_vessel_data.count(key) == 1)
                {
                    pVesselInfoVector[idx]->SetValue(vessel_index, vtk_vessel_data[key]);
                }
    //            // Otherwise check the generic data
    //            else if(generic_vessel_data.count(key) == 1)
    //            {
    //                if(generic_vessel_data[key].type() == typeid(double))
    //                {
    //                    double cast_value = boost::any_cast<double>(generic_vessel_data[key]);
    //                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
    //                }
    //                else if(generic_vessel_data[key].type() == typeid(unsigned))
    //                {
    //                    double cast_value = double(boost::any_cast<unsigned>(generic_vessel_data[key]));
    //                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
    //                }
    //                else if(generic_vessel_data[key].type() == typeid(bool))
    //                {
    //                    double cast_value = double(boost::any_cast<bool>(generic_vessel_data[key]));
    //                    pVesselInfoVector[idx]->SetValue(vessel_index, cast_value);
    //                }
    //            }
            }
            vessel_index++;
        }
        pPolyData->SetPoints(pPoints);
        pPolyData->SetLines(pLines);

        for (unsigned i = 0; i < pVesselInfoVector.size(); i++)
        {
            pPolyData->GetCellData()->AddArray(pVesselInfoVector[i]);
        }

        for (unsigned i = 0; i < pNodeInfoVector.size(); i++)
        {
            pPolyData->GetPointData()->AddArray(pNodeInfoVector[i]);
        }
    }

    return pPolyData;
}

template <unsigned DIM>
void VtkVesselNetworkWriter<DIM>::Write()
{
    if(mFilename.empty())
    {
        EXCEPTION("No file name set for VtkVesselNetworkWriter");
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mFilename.c_str());
    writer->SetInput(this->GetOutput());
    writer->Write();
}

// Explicit instantiation
template class VesselNetwork<2>;
template class VesselNetwork<3>;

