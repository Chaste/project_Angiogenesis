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
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include "vtkXMLPolyDataWriter.h"
#include "vtkCleanPolyData.h"
#include "vtkSelectEnclosedPoints.h"
#include "vtkTriangleFilter.h"
#include "vtkSTLWriter.h"
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "VesselSurfaceGenerator.hpp"

#include "Part.hpp"

struct triangulateio;

template<unsigned DIM>
Part<DIM>::Part() :
        mFacets(),
        mVtkPart(vtkSmartPointer<vtkPolyData>()),
        mHoleMarkers(),
        mRegionMarkers()
{
}

template<unsigned DIM>
boost::shared_ptr<Part<DIM> > Part<DIM>::Create()
{
    MAKE_PTR(Part<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
Part<DIM>::~Part()
{

}

template<unsigned DIM>
boost::shared_ptr<Polygon> Part<DIM>::AddCircle(double radius, c_vector<double, DIM> centre, unsigned numSegments)
{
    std::vector<boost::shared_ptr<Vertex> > vertices;
    double seg_angle = 2.0 * M_PI / double(numSegments);
    for (unsigned idx = 0; idx < numSegments; idx++)
    {
        double angle = seg_angle * double(idx);
        double x = radius * std::cos(angle) + centre[0];
        double y = radius * std::sin(angle) + centre[1];
        if(DIM==3)
        {
            vertices.push_back(Vertex::Create(x, y, centre[2]));
        }
        else
        {
            vertices.push_back(Vertex::Create(x, y, 0.0));
        }
    }
    return AddPolygon(vertices);
}

template<unsigned DIM>
void Part<DIM>::AddCylinder(double radius, double depth, c_vector<double, DIM> centre, unsigned numSegments)
{
    boost::shared_ptr<Polygon> p_circle = AddCircle(radius, centre, numSegments);
    Extrude(p_circle, depth);
}

template<unsigned DIM>
void Part<DIM>::AddCuboid(double sizeX, double sizeY, double sizeZ, c_vector<double, DIM> origin)
{
    boost::shared_ptr<Polygon> p_rectangle = AddRectangle(sizeX, sizeY, origin);
    Extrude(p_rectangle, sizeZ);
}

template<unsigned DIM>
void Part<DIM>::AddHoleMarker(c_vector<double, DIM> hole)
{
    mHoleMarkers.push_back(hole);
}

template<unsigned DIM>
boost::shared_ptr<Polygon> Part<DIM>::AddPolygon(std::vector<boost::shared_ptr<Vertex> > vertices, bool newFacet,
                                                                                   boost::shared_ptr<Facet> pFacet)
{
    boost::shared_ptr<Polygon> p_polygon = Polygon::Create(vertices);
    AddPolygon(p_polygon, newFacet, pFacet);
    return p_polygon;
}

template<unsigned DIM>
boost::shared_ptr<Polygon> Part<DIM>::AddPolygon(boost::shared_ptr<Polygon> pPolygon, bool newFacet, boost::shared_ptr<Facet> pFacet)
{
    if (!pFacet)
    {
        if (mFacets.size() == 0 || newFacet)
        {
            mFacets.push_back(Facet::Create(pPolygon));
        }
        else
        {
            mFacets[0]->AddPolygon(pPolygon);
        }
    }
    else
    {
        pFacet->AddPolygon(pPolygon);
    }

    return pPolygon;
}

template<unsigned DIM>
boost::shared_ptr<Polygon> Part<DIM>::AddRectangle(double sizeX, double sizeY, c_vector<double, DIM> origin)
{
    std::vector<boost::shared_ptr<Vertex> > vertices;
    if(DIM==3)
    {
        vertices.push_back(Vertex::Create(origin[0], origin[1], origin[2]));
        vertices.push_back(Vertex::Create(origin[0] + sizeX, origin[1], origin[2]));
        vertices.push_back(Vertex::Create(origin[0] + sizeX, origin[1] + sizeY, origin[2]));
        vertices.push_back(Vertex::Create(origin[0], origin[1] + sizeY, origin[2]));
    }
    else
    {
        vertices.push_back(Vertex::Create(origin[0], origin[1], 0.0));
        vertices.push_back(Vertex::Create(origin[0] + sizeX, origin[1], 0.0));
        vertices.push_back(Vertex::Create(origin[0] + sizeX, origin[1] + sizeY, 0.0));
        vertices.push_back(Vertex::Create(origin[0], origin[1] + sizeY, 0.0));
    }
    return AddPolygon(vertices);
}

template<unsigned DIM>
void Part<DIM>::AddVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pVesselNetwork, bool surface)
{
    if (!surface)
    {
        std::vector<boost::shared_ptr<Vertex> > vertices;
        std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = pVesselNetwork->GetNodes();
        for (unsigned idx = 0; idx < nodes.size(); idx++)
        {
            vertices.push_back(Vertex::Create(nodes[idx]->GetLocation().rGetLocation()));
        }

        // If vertices lie on any existing facets add the vertex to the facet
        for (unsigned kdx = 0; kdx < vertices.size(); kdx++)
        {
            for (unsigned idx = 0; idx < mFacets.size(); idx++)
            {
                if(mFacets[idx]->ContainsPoint(vertices[kdx]->rGetLocation()))
                {
                    mFacets[idx]->AddPolygon(Polygon::Create(vertices[kdx]));
                }
            }
        }

        // Create polygons and facets for each vessel
        std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = pVesselNetwork->GetVessels();
        for (unsigned idx = 0; idx < vessels.size(); idx++)
        {
            std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vessels[idx]->GetSegments();
            std::vector<boost::shared_ptr<Vertex> > segment_vertices;
            for (unsigned jdx = 0; jdx < segments.size(); jdx++)
            {
                unsigned node0_index = pVesselNetwork->GetNodeIndex(segments[jdx]->GetNode(0));
                unsigned node1_index = pVesselNetwork->GetNodeIndex(segments[jdx]->GetNode(1));
                if (jdx == 0)
                {
                    segment_vertices.push_back(vertices[node0_index]);
                }
                segment_vertices.push_back(vertices[node1_index]);
            }
            boost::shared_ptr<Polygon> p_polygon = Polygon::Create(segment_vertices);
            mFacets.push_back(Facet::Create(p_polygon));
        }
    }
    else
    {
        // Add any polygons on existing facets to the facet
        VesselSurfaceGenerator<DIM> generator(pVesselNetwork);
        std::vector<boost::shared_ptr<Polygon> > polygons = generator.GetSurfacePolygons();
        std::vector<bool> polygon_on_facet;

        for (unsigned idx = 0; idx < polygons.size(); idx++)
        {
            bool on_facet = false;
            c_vector<double, DIM> poly_centroid = polygons[idx]->GetCentroid();
            for (unsigned jdx = 0; jdx < mFacets.size(); jdx++)
            {
                if (mFacets[jdx]->ContainsPoint(poly_centroid))
                {
                    on_facet = true;
                    mFacets[jdx]->AddPolygon(polygons[idx]);
                }
            }
            polygon_on_facet.push_back(on_facet);
        }

        // Create polygons and facets for each vessel
        for (unsigned idx = 0; idx < polygons.size(); idx++)
        {
            if (!polygon_on_facet[idx])
            {
                mFacets.push_back(Facet::Create(polygons[idx]));
            }
        }

        std::vector<c_vector<double, DIM> > hole_locations = generator.GetHoles();
        for(unsigned idx=0; idx<hole_locations.size(); idx++)
        {
            AddHoleMarker(hole_locations[idx]);
        }
    }

}

template<unsigned DIM>
void Part<DIM>::Extrude(boost::shared_ptr<Polygon> pPolygon, double depth)
{
    if(DIM==2)
    {
        EXCEPTION("Only parts in 3D space can be extruded.");
    }
    // Loop through the vertices and create new ones at the offset depth
    std::vector<boost::shared_ptr<Vertex> > original_vertices = pPolygon->GetVertices();
    std::vector<boost::shared_ptr<Vertex> > new_vertices;
    for (unsigned idx = 0; idx < original_vertices.size(); idx++)
    {
        c_vector<double, DIM> location = original_vertices[idx]->rGetLocation();
        new_vertices.push_back(Vertex::Create(location[0], location[1], location[2] + depth));
    }

    // Every straight edge is now a planar face, with 3 new edges ordered in CCW
    for (unsigned idx = 0; idx < original_vertices.size(); idx++)
    {
        unsigned index2;
        if (idx != original_vertices.size() - 1)
        {
            index2 = idx + 1;
        }
        else
        {
            index2 = 0;
        }
        boost::shared_ptr<Polygon> p_polygon = Polygon::Create(original_vertices[idx]);
        p_polygon->AddVertex(original_vertices[index2]);
        p_polygon->AddVertex(new_vertices[index2]);
        p_polygon->AddVertex(new_vertices[idx]);
        mFacets.push_back(Facet::Create(p_polygon));
    }

    // Close the lid
    boost::shared_ptr<Polygon> p_polygon = Polygon::Create(new_vertices);
    mFacets.push_back(Facet::Create(p_polygon));
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > Part<DIM>::GetHoleMarkers()
{
    return mHoleMarkers;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<Vertex> > Part<DIM>::GetVertices()
{
    std::vector<boost::shared_ptr<Polygon> > polygons = GetPolygons();
    std::set<boost::shared_ptr<Vertex> > unique_vertices;
    for (unsigned idx = 0; idx < polygons.size(); idx++)
    {
        std::vector<boost::shared_ptr<Vertex> > polygon_vertices = polygons[idx]->GetVertices();
        std::copy(polygon_vertices.begin(), polygon_vertices.end(),
                  std::inserter(unique_vertices, unique_vertices.end()));
    }

    std::vector<boost::shared_ptr<Vertex> > vertices;
    vertices.insert(vertices.end(), unique_vertices.begin(), unique_vertices.end());

    for (unsigned idx = 0; idx < vertices.size(); idx++)
    {
        vertices[idx]->SetIndex(idx);
    }
    return vertices;
}

template<unsigned DIM>
std::vector<unsigned> Part<DIM>::GetContainingGridIndices(unsigned num_x, unsigned num_y, unsigned num_z, double spacing)
{
    std::vector<unsigned> location_indices;
    for(unsigned kdx=0; kdx<num_z; kdx++)
    {
        for(unsigned jdx=0; jdx<num_y; jdx++)
        {
            for(unsigned idx=0; idx<num_x; idx++)
            {
                c_vector<double,3> location;
                location[0] = double(idx) * spacing;
                location[1] = double(jdx) * spacing;
                location[2] = double(kdx) * spacing;
                unsigned index = idx + num_x * jdx + num_x * num_y * kdx;
                bool update = false;
                if(index==0)
                {
                    update = true;
                }
                if(IsPointInPart(location, update))
                {
                    location_indices.push_back(index);
                }
            }
        }
    }
    return location_indices;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > Part<DIM>::GetVertexLocations()
{
    std::vector<boost::shared_ptr<Vertex> > vertices = GetVertices();
    std::vector<c_vector<double, DIM> > locations;

    typename std::vector<boost::shared_ptr<Vertex> >::iterator iter;
    for (iter = vertices.begin(); iter != vertices.end(); iter++)
    {
        c_vector<double, DIM> this_location;
        for(unsigned jdx=0; jdx<DIM;jdx++)
        {
            this_location[jdx] = (*iter)->rGetLocation()[jdx];
        }
        locations.push_back(this_location);
    }
    return locations;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<Polygon> > Part<DIM>::GetPolygons()
{
    std::vector<boost::shared_ptr<Polygon> > polygons;
    for (unsigned idx = 0; idx < mFacets.size(); idx++)
    {
        std::vector<boost::shared_ptr<Polygon> > facet_polygons = mFacets[idx]->GetPolygons();
        polygons.insert(polygons.end(), facet_polygons.begin(), facet_polygons.end());
    }
    return polygons;
}

template<unsigned DIM>
c_vector<double, 2*DIM> Part<DIM>::GetBoundingBox()
{
    std::vector<boost::shared_ptr<Vertex> > vertices = GetVertices();
    c_vector<double, 2*DIM> box;

    for (unsigned idx = 0; idx < vertices.size(); idx++)
    {
        for (unsigned jdx = 0; jdx < DIM; jdx++)
        {
            if (idx == 0)
            {
                box[2 * jdx] = vertices[idx]->rGetLocation()[jdx];
                box[2 * jdx + 1] = vertices[idx]->rGetLocation()[jdx];
            }
            else
            {
                if (vertices[idx]->rGetLocation()[jdx] < box[2 * jdx])
                {
                    box[2 * jdx] = vertices[idx]->rGetLocation()[jdx];
                }
                if (vertices[idx]->rGetLocation()[jdx] > box[2 * jdx + 1])
                {
                    box[2 * jdx + 1] = vertices[idx]->rGetLocation()[jdx];
                }
            }
        }
    }
    return box;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<Facet> > Part<DIM>::GetFacets()
{
    return mFacets;
}

template<unsigned DIM>
std::vector<std::pair<unsigned, unsigned> > Part<DIM>::GetSegmentIndices()
{
    // Make sure the vertex indexes are up-to-date.
    GetVertices();

    std::vector<std::pair<unsigned, unsigned> > indexes;
    std::vector<boost::shared_ptr<Polygon> > polygons = mFacets[0]->GetPolygons();

    for (unsigned idx = 0; idx < polygons.size(); idx++)
    {
        if (polygons[idx]->GetVertices().size() == 2)
        {
            indexes.push_back(
                    std::pair<unsigned, unsigned>(polygons[idx]->GetVertices()[0]->GetIndex(),
                                                  polygons[idx]->GetVertices()[1]->GetIndex()));
        }
        else if (polygons[idx]->GetVertices().size() > 1)
        {
            std::vector<boost::shared_ptr<Vertex> > vertices = polygons[idx]->GetVertices();
            for (unsigned jdx = 0; jdx < vertices.size() - 1; jdx++)
            {
                indexes.push_back(
                        std::pair<unsigned, unsigned>(vertices[jdx]->GetIndex(), vertices[jdx + 1]->GetIndex()));
            }
            indexes.push_back(
                    std::pair<unsigned, unsigned>(vertices[vertices.size() - 1]->GetIndex(), vertices[0]->GetIndex()));

        }
    }
    return indexes;
}

template<unsigned DIM>
vtkSmartPointer<vtkPolyData> Part<DIM>::GetVtk(bool update)
{
    vtkSmartPointer<vtkPolyData> p_part_data = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> p_vertices = vtkSmartPointer<vtkPoints>::New();

    p_part_data->Allocate(1, 1);
    // Loop through each polygon, collect the vertices and set correct point ids
    std::vector<boost::shared_ptr<Polygon> > polygons = GetPolygons();
    unsigned vert_counter = 0;
    for (vtkIdType idx = 0; idx < vtkIdType(polygons.size()); idx++)
    {
        std::vector<boost::shared_ptr<Vertex> > vertices = polygons[idx]->GetVertices();
        vtkSmartPointer<vtkPolygon> p_polygon = vtkSmartPointer<vtkPolygon>::New();
        p_polygon->GetPointIds()->SetNumberOfIds(vertices.size());
        for (vtkIdType jdx = 0; jdx < vtkIdType(vertices.size()); jdx++)
        {
            c_vector<double, 3> location = vertices[jdx]->rGetLocation();
            if(DIM==3)
            {
                p_vertices->InsertNextPoint(location[0], location[1], location[2]);
            }
            else
            {
                p_vertices->InsertNextPoint(location[0], location[1], 0.0);
            }

            p_polygon->GetPointIds()->SetId(jdx, vert_counter);
            vert_counter++;
        }
        p_part_data->InsertNextCell(p_polygon->GetCellType(), p_polygon->GetPointIds());
    }
    p_part_data->SetPoints(p_vertices);
    vtkSmartPointer<vtkCleanPolyData> p_clean_data = vtkSmartPointer<vtkCleanPolyData>::New();
    p_clean_data->SetInput(p_part_data);
    mVtkPart = p_clean_data->GetOutput();
    return mVtkPart;
}

template<unsigned DIM>
bool Part<DIM>::IsPointInPart(c_vector<double, DIM> location, bool update)
{
    vtkSmartPointer<vtkPolyData> p_part = GetVtk(update);
    vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();

    if(DIM==3)
    {
        p_points->InsertNextPoint(location[0], location[1], location[2]);
    }
    else
    {
        p_points->InsertNextPoint(location[0], location[1], 0.0);
    }

    vtkSmartPointer<vtkPolyData> p_point_data = vtkSmartPointer<vtkPolyData>::New();
    p_point_data->SetPoints(p_points);

    //Points inside test
    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    #if VTK_MAJOR_VERSION <= 5
    selectEnclosedPoints->SetInput(p_point_data);
    #else
    selectEnclosedPoints->SetInputData(p_point_data);
    #endif
    #if VTK_MAJOR_VERSION <= 5
    selectEnclosedPoints->SetSurface(p_part);
    #else
    selectEnclosedPoints->SetSurfaceData(p_part);
    #endif
    selectEnclosedPoints->Update();

    return selectEnclosedPoints->IsInside(0);
}

template<unsigned DIM>
void Part<DIM>::MergeCoincidentVertices()
{
    // Loop through the nodes of each polygon. If it is in another polygon, replace it.
    std::vector<boost::shared_ptr<Polygon> > polygons = GetPolygons();
    for(unsigned idx=0; idx<polygons.size(); idx++)
    {
        for(unsigned jdx=0; jdx<polygons.size(); jdx++)
        {
            if(idx != jdx)
            {
                std::vector<boost::shared_ptr<Vertex> > p1_verts = polygons[idx]->GetVertices();
                std::vector<boost::shared_ptr<Vertex> > p2_verts = polygons[jdx]->GetVertices();
                for(unsigned mdx=0; mdx<p1_verts.size(); mdx++)
                {
                    for(unsigned ndx=0; ndx<p2_verts.size(); ndx++)
                    {
                        if(norm_2(p2_verts[ndx]->rGetLocation()- p1_verts[mdx]->rGetLocation())< 1.e-6)
                        {
                            polygons[jdx]->ReplaceVertex(ndx, polygons[idx]->GetVertex(mdx));
                        }
                    }
                }
            }
        }
    }

    for(unsigned idx=0; idx<mFacets.size(); idx++)
    {
        mFacets[idx]->UpdateVertices();
    }
}

template<unsigned DIM>
void Part<DIM>::Translate(c_vector<double, DIM> vector)
{
    std::vector<boost::shared_ptr<Vertex> > vertices = GetVertices();
    {
        for (unsigned idx = 0; idx < vertices.size(); idx++)
        {
            vertices[idx]->Translate(vector);
        }
    }
}

template<unsigned DIM>
void Part<DIM>::Write(const std::string& fileName)
{
    GetVtk();
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(mVtkPart);
    writer->Write();
}

template<unsigned DIM>
void Part<DIM>::WriteStl(const std::string& rFilename)
{
    vtkSmartPointer<vtkTriangleFilter> p_tri_filter = vtkSmartPointer<vtkTriangleFilter>::New();
    p_tri_filter->SetInput(GetVtk());

    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName(rFilename.c_str());
    stlWriter->SetInput(p_tri_filter->GetOutput());
    stlWriter->SetFileTypeToASCII();
    stlWriter->Write();
}

// Explicit instantiation
template class Part<2> ;
template class Part<3> ;
