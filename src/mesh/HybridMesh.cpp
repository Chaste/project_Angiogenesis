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

#include <boost/lexical_cast.hpp>
#include "Exception.hpp"
#include "Warnings.hpp"
#include "Facet.hpp"
#include "Polygon.hpp"
#include "HybridMesh.hpp"
#include "Element.hpp"
#include "UblasVectorInclude.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HybridMesh<ELEMENT_DIM, SPACE_DIM>::HybridMesh()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<HybridMesh<ELEMENT_DIM, SPACE_DIM> > HybridMesh<ELEMENT_DIM, SPACE_DIM>::Create()
{
    MAKE_PTR(HybridMesh<ELEMENT_DIM>, pSelf);
    return pSelf;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HybridMesh<ELEMENT_DIM, SPACE_DIM>::~HybridMesh()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::Mesh2d(boost::shared_ptr<Part<SPACE_DIM> > pPart, double maxElementArea)
{
    std::vector<c_vector<double, SPACE_DIM> > vertex_locations = pPart->GetVertexLocations();
    unsigned num_vertices = vertex_locations.size();
    struct triangulateio mesher_input, mesher_output;
    this->InitialiseTriangulateIo(mesher_input);
    this->InitialiseTriangulateIo(mesher_output);

    mesher_input.pointlist = (double *) malloc(num_vertices * 2 * sizeof(double));
    mesher_input.numberofpoints = num_vertices;
    for (unsigned idx = 0; idx < num_vertices; idx++)
    {
        for (unsigned jdx = 0; jdx < 2; jdx++)
        {
            mesher_input.pointlist[2 * idx + jdx] = vertex_locations[idx][jdx];
        }
    }
    std::vector<std::pair<unsigned, unsigned> > segments = pPart->GetSegmentIndices();
    unsigned num_segments = segments.size();

    mesher_input.segmentlist = (int *) malloc(num_segments * 2 * sizeof(int));
    mesher_input.numberofsegments = num_segments;
    for (unsigned idx = 0; idx < num_segments; idx++)
    {
        mesher_input.segmentlist[2 * idx] = int(segments[idx].first);
        mesher_input.segmentlist[2 * idx + 1] = int(segments[idx].second);
    }

    std::string mesher_command = "pqQze";
    if (maxElementArea > 0.0)
    {
        mesher_command += "a" + boost::lexical_cast<std::string>(maxElementArea);
    }
    triangulate((char*) mesher_command.c_str(), &mesher_input, &mesher_output, NULL);

    this->ImportFromMesher(mesher_output, mesher_output.numberoftriangles, mesher_output.trianglelist,
                           mesher_output.numberofedges, mesher_output.edgelist, mesher_output.edgemarkerlist);

    //Tidy up triangle
    this->FreeTriangulateIo(mesher_input);
    this->FreeTriangulateIo(mesher_output);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::Mesh3d(boost::shared_ptr<Part<SPACE_DIM> > pPart, double maxElementArea)
{
    std::vector<c_vector<double, SPACE_DIM> > vertex_locations = pPart->GetVertexLocations();
    std::vector<c_vector<double, SPACE_DIM> > hole_locations = pPart->GetHoleMarkers();
    unsigned num_vertices = vertex_locations.size();
    unsigned num_holes = hole_locations.size();
    std::vector<boost::shared_ptr<Facet> > facets = pPart->GetFacets();
    unsigned num_facets = facets.size();

    class tetgen::tetgenio mesher_input, mesher_output;

    tetgen::tetgenio::facet *f;
    tetgen::tetgenio::polygon *p;
    mesher_input.pointlist = new double[(num_vertices) * 3];
    mesher_input.numberofpoints = num_vertices;

    for (unsigned idx = 0; idx < num_vertices; idx++)
    {
        for (unsigned jdx = 0; jdx < 3; jdx++)
        {
            mesher_input.pointlist[3 * idx + jdx] = vertex_locations[idx][jdx];
        }
    }

    // Add the holes
    mesher_input.holelist = new double[(num_holes) * 3];
    mesher_input.numberofholes = num_holes;
    for (unsigned idx = 0; idx < num_holes; idx++)
    {
        for (unsigned jdx = 0; jdx < 3; jdx++)
        {
            mesher_input.holelist[3 * idx + jdx] = hole_locations[idx][jdx];
        }
    }

    mesher_input.numberoffacets = num_facets;
    mesher_input.facetlist = new tetgen::tetgenio::facet[num_facets];
    mesher_input.facetmarkerlist = new int[num_facets];
    for (unsigned idx = 0; idx < num_facets; idx++)
    {
        mesher_input.facetmarkerlist[idx] = 0;
        f = &mesher_input.facetlist[idx];
        std::vector<boost::shared_ptr<Polygon> > polygons = facets[idx]->GetPolygons();
        f->numberofpolygons = polygons.size();
        f->polygonlist = new tetgen::tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        for (unsigned jdx = 0; jdx < polygons.size(); jdx++)
        {
            p = &f->polygonlist[jdx];
            p->numberofvertices = polygons[jdx]->GetVertices().size();
            p->vertexlist = new int[p->numberofvertices];
            for (unsigned kdx = 0; kdx < polygons[jdx]->GetVertices().size(); kdx++)
            {
                p->vertexlist[kdx] = int(polygons[jdx]->GetVertices()[kdx]->GetIndex());
            }
        }
    }
    std::string mesher_command = "pqQz";
    if (maxElementArea > 0.0)
    {
        mesher_command += "a" + boost::lexical_cast<std::string>(maxElementArea);
    }

    // Library call
    tetgen::tetrahedralize((char*) mesher_command.c_str(), &mesher_input, &mesher_output);

    this->ImportFromTetgen(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist,
                           mesher_output.numberoftrifaces, mesher_output.trifacelist, NULL);

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::GenerateFromStl(const std::string& filename, double maxElementArea,
                                                         std::vector<c_vector<double, SPACE_DIM> > holes)
{
    class tetgen::tetgenio mesher_input, mesher_output;
    char * writable = new char[filename.size() + 1];
    std::copy(filename.begin(), filename.end(), writable);
    writable[filename.size()] = '\0';
    mesher_input.load_stl(writable);

    unsigned num_holes = holes.size();
    mesher_input.holelist = new double[(num_holes) * 3];
    mesher_input.numberofholes = num_holes;
    for (unsigned idx = 0; idx < num_holes; idx++)
    {
        for (unsigned jdx = 0; jdx < 3; jdx++)
        {
            mesher_input.holelist[3 * idx + jdx] = holes[idx][jdx];
        }
    }

    std::string mesher_command = "pqQz";
    if (maxElementArea > 0.0)
    {
        mesher_command += "a" + boost::lexical_cast<std::string>(maxElementArea);
    }

    // Library call
    tetgen::tetrahedralize((char*) mesher_command.c_str(), &mesher_input, &mesher_output);

    this->ImportFromTetgen(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist,
                           mesher_output.numberoftrifaces, mesher_output.trifacelist, NULL);

    delete[] writable;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::ImportFromTetgen(tetgen::tetgenio& mesherOutput, unsigned numberOfElements,
                                                          int *elementList, unsigned numberOfFaces, int *faceList,
                                                          int *edgeMarkerList)
{
    unsigned nodes_per_element = mesherOutput.numberofcorners;

    assert(nodes_per_element == ELEMENT_DIM + 1 || nodes_per_element == (ELEMENT_DIM + 1) * (ELEMENT_DIM + 2) / 2);

    for (unsigned i = 0; i < this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryElements.clear();
    this->mBoundaryNodes.clear();

    // Construct the nodes
    for (unsigned node_index = 0; node_index < (unsigned) mesherOutput.numberofpoints; node_index++)
    {
        this->mNodes.push_back(new Node<SPACE_DIM>(node_index, &mesherOutput.pointlist[node_index * SPACE_DIM], false));
    }

    // Construct the elements
    this->mElements.reserve(numberOfElements);

    unsigned real_element_index = 0;
    for (unsigned element_index = 0; element_index < numberOfElements; element_index++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned j = 0; j < ELEMENT_DIM + 1; j++)
        {
            unsigned global_node_index = elementList[element_index * (nodes_per_element) + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);

        }

        /*
         * For some reason, tetgen in library mode makes its initial Delaunay mesh
         * with very thin slivers. Hence we expect to ignore some of the elements!
         */
        Element<ELEMENT_DIM, SPACE_DIM>* p_element;
        try
        {
            p_element = new Element<ELEMENT_DIM, SPACE_DIM>(real_element_index, nodes);

            // Shouldn't throw after this point
            this->mElements.push_back(p_element);

            // Add the internals to quadratics
            for (unsigned j = ELEMENT_DIM + 1; j < nodes_per_element; j++)
            {
                unsigned global_node_index = elementList[element_index * nodes_per_element + j];
                assert(global_node_index < this->mNodes.size());
                this->mElements[real_element_index]->AddNode(this->mNodes[global_node_index]);
                this->mNodes[global_node_index]->AddElement(real_element_index);
                this->mNodes[global_node_index]->MarkAsInternal();
            }
            real_element_index++;
        }
        catch (Exception &)
        {
            if (SPACE_DIM == 2)
            {
                WARNING("Triangle has produced a zero area (collinear) element");
            }
            else
            {
                WARNING("Tetgen has produced a zero volume (coplanar) element");
            }
        }
    }

    // Construct the BoundaryElements (and mark boundary nodes)
    unsigned next_boundary_element_index = 0;
    for (unsigned boundary_element_index = 0; boundary_element_index < numberOfFaces; boundary_element_index++)
    {
        /*
         * Tetgen produces only boundary faces (set edgeMarkerList to NULL).
         * Triangle marks which edges are on the boundary.
         */
        if (edgeMarkerList == NULL || edgeMarkerList[boundary_element_index] == 1)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j = 0; j < ELEMENT_DIM; j++)
            {
                unsigned global_node_index = faceList[boundary_element_index * ELEMENT_DIM + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
            }

            /*
             * For some reason, tetgen in library mode makes its initial Delaunay mesh
             * with very thin slivers. Hence we expect to ignore some of the elements!
             */
            BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_b_element;
            try
            {
                p_b_element = new BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>(next_boundary_element_index, nodes);
                this->mBoundaryElements.push_back(p_b_element);
                next_boundary_element_index++;
            }
            catch (Exception &)
            {
                // Tetgen is feeding us lies  //Watch this space for coverage
                assert(SPACE_DIM == 3);
            }
        }
    }

    this->RefreshJacobianCachedData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::GenerateFromPart(boost::shared_ptr<Part<SPACE_DIM> > pPart,
                                                          double maxElementArea)
{
    // For 2D parts use triangle
    if (ELEMENT_DIM == 2)
    {
        c_vector<double, 2 * ELEMENT_DIM> bounding_box = pPart->GetBoundingBox();
        if (SPACE_DIM == 2)
        {
            Mesh2d(pPart, maxElementArea);
        }
        else if (std::abs(bounding_box[4]) < 1.e-6 && std::abs(bounding_box[5]) < 1.e-6)
        {
            Mesh2d(pPart, maxElementArea);
        }
        else
        {
            EXCEPTION("For now 2D meshing is only supported for parts with z=0.");
        }
    }
    // Try to use tetgen
    else
    {
        c_vector<double, 2 * ELEMENT_DIM> bounding_box = pPart->GetBoundingBox();
        if (std::abs(bounding_box[4]) < 1.e-6 && std::abs(bounding_box[5]) < 1.e-6)
        {
            EXCEPTION("The part is two-dimensional, use the 2D meshing functionality.");
        }
        else
        {
            Mesh3d(pPart, maxElementArea);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > HybridMesh<ELEMENT_DIM, SPACE_DIM>::GetConnectivity()
{
    std::vector<std::vector<unsigned> > connectivity;
    unsigned num_elements = AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements();
    for (unsigned idx = 0; idx < num_elements; idx++)
    {
        std::vector<unsigned> node_indexes;
        unsigned num_nodes = AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(idx)->GetNumNodes();
        for (unsigned jdx = 0; jdx < num_nodes; jdx++)
        {
            node_indexes.push_back(
                    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(idx)->GetNodeGlobalIndex(jdx));
        }
        connectivity.push_back(node_indexes);
    }
    return connectivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<double> > HybridMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeLocations()
{
    std::vector<std::vector<double> > locations;
    unsigned num_nodes = AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes();
    for (unsigned idx = 0; idx < num_nodes; idx++)
    {
        c_vector<double, SPACE_DIM> location =
                AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(idx)->rGetLocation();
        std::vector<double> vec_location;
        for (unsigned jdx = 0; jdx < SPACE_DIM; jdx++)
        {
            vec_location.push_back(location[jdx]);
        }
        locations.push_back(vec_location);
    }
    return locations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::InitialiseTriangulateIo(triangulateio& mesherIo)
{
    mesherIo.numberofpoints = 0;
    mesherIo.pointlist = NULL;
    mesherIo.numberofpointattributes = 0;
    mesherIo.pointattributelist = (double *) NULL;
    mesherIo.pointmarkerlist = (int *) NULL;
    mesherIo.segmentlist = NULL;
    mesherIo.segmentmarkerlist = (int *) NULL;
    mesherIo.numberofsegments = 0;
    mesherIo.numberofholes = 0;
    mesherIo.numberofregions = 0;
    mesherIo.trianglelist = (int *) NULL;
    mesherIo.triangleattributelist = (double *) NULL;
    mesherIo.numberoftriangleattributes = 0;
    mesherIo.edgelist = (int *) NULL;
    mesherIo.edgemarkerlist = (int *) NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HybridMesh<ELEMENT_DIM, SPACE_DIM>::FreeTriangulateIo(triangulateio& mesherIo)
{
    if (mesherIo.numberofpoints != 0)
    {
        mesherIo.numberofpoints = 0;
        free(mesherIo.pointlist);
    }

    if (mesherIo.numberofsegments != 0)
    {
        mesherIo.numberofsegments = 0;
        free(mesherIo.segmentlist);
    }

    // These (and the above) should actually be safe since we explicity set to NULL above
    free(mesherIo.pointattributelist);
    free(mesherIo.pointmarkerlist);
    free(mesherIo.segmentmarkerlist);
    free(mesherIo.trianglelist);
    free(mesherIo.triangleattributelist);
    free(mesherIo.edgelist);
    free(mesherIo.edgemarkerlist);
}

// Explicit instantiation
template class HybridMesh<2> ;
template class HybridMesh<3> ;
