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
#include "Polygon.hpp"
#include "Vertex.hpp"

#include "MappableGridGenerator.hpp"

template<unsigned DIM>
MappableGridGenerator<DIM>::MappableGridGenerator()
{
}

template<unsigned DIM>
MappableGridGenerator<DIM>::~MappableGridGenerator()
{

}

template<unsigned DIM>
boost::shared_ptr<Part<DIM> > MappableGridGenerator<DIM>::GeneratePlane(unsigned numX, unsigned numY)
{
    // Make a regular grid of polygons
    // Make the vertices
    std::vector<boost::shared_ptr<Vertex> > vertices;
    for(unsigned jdx=0; jdx< numY; jdx++)
    {
        for(unsigned idx=0; idx<numX; idx++)
        {
            double x_location = double(idx);
            double y_location = double(jdx);
            double z_location = 0.0;

            // Create the vertices
            vertices.push_back(Vertex::Create(x_location, y_location, z_location));

        }
    }

    for(unsigned jdx=0; jdx< numY; jdx++)
    {
        for(unsigned idx=0; idx<numX; idx++)
        {
            double x_location = double(idx);
            double y_location = double(jdx);
            double z_location = 1.0;

            // Create the vertices
            vertices.push_back(Vertex::Create(x_location, y_location, z_location));

        }
    }

    // Make the polygons
    // Front face
    std::vector<boost::shared_ptr<Polygon> > polygons;
    for(unsigned jdx=0; jdx< numY - 1; jdx++)
    {
        for(unsigned idx=0; idx<numX - 1; idx++)
        {
            unsigned front_left_index = idx + numX * jdx;
            unsigned front_right_index = idx + 1 + numX * jdx;
            unsigned front_left_top_index = idx + numX * (jdx+1);
            unsigned front_right_top_index = idx + 1 + numX * (jdx+1);

            std::vector<boost::shared_ptr<Vertex> > poly_vertices;
            poly_vertices.push_back(vertices[front_left_index]);
            poly_vertices.push_back(vertices[front_right_index]);
            poly_vertices.push_back(vertices[front_right_top_index]);
            poly_vertices.push_back(vertices[front_left_top_index]);
            polygons.push_back(Polygon::Create(poly_vertices));
        }
    }

    // Back face
    for(unsigned jdx=0; jdx< numY - 1; jdx++)
    {
        for(unsigned idx=0; idx<numX - 1; idx++)
        {
            unsigned front_left_index = idx + numX * jdx + numX*numY;
            unsigned front_right_index = idx + 1 + numX * jdx + numX*numY;
            unsigned front_left_top_index = idx + numX * (jdx+1) + numX*numY;
            unsigned front_right_top_index = idx + 1 + numX * (jdx+1) + numX*numY;

            std::vector<boost::shared_ptr<Vertex> > poly_vertices;
            poly_vertices.push_back(vertices[front_left_index]);
            poly_vertices.push_back(vertices[front_right_index]);
            poly_vertices.push_back(vertices[front_right_top_index]);
            poly_vertices.push_back(vertices[front_left_top_index]);
            polygons.push_back(Polygon::Create(poly_vertices));
        }
    }

    // Left face
    for(unsigned jdx=0; jdx< numY - 1; jdx++)
    {
        unsigned front_index = numX * jdx;
        unsigned top_front_index = numX * (jdx+1);
        unsigned back_index = numX * jdx + numX*numY;
        unsigned top_back_index = numX * (jdx+1) + numX*numY;

        std::vector<boost::shared_ptr<Vertex> > poly_vertices;
        poly_vertices.push_back(vertices[front_index]);
        poly_vertices.push_back(vertices[top_front_index]);
        poly_vertices.push_back(vertices[top_back_index]);
        poly_vertices.push_back(vertices[back_index]);
        polygons.push_back(Polygon::Create(poly_vertices));
    }

    // Right face
    for(unsigned jdx=0; jdx< numY - 1; jdx++)
    {
        unsigned front_index = numX * (jdx+1) - 1;
        unsigned top_front_index = numX * (jdx+2) - 1;
        unsigned back_index = numX * (jdx + 1) - 1 + numX*numY;
        unsigned top_back_index = numX * (jdx+2) -1 + numX*numY;

        std::vector<boost::shared_ptr<Vertex> > poly_vertices;
        poly_vertices.push_back(vertices[front_index]);
        poly_vertices.push_back(vertices[top_front_index]);
        poly_vertices.push_back(vertices[top_back_index]);
        poly_vertices.push_back(vertices[back_index]);
        polygons.push_back(Polygon::Create(poly_vertices));
    }

    // Bottom face
    for(unsigned idx=0; idx< numX - 1; idx++)
    {
        unsigned front_index = idx;
        unsigned front_right_index = idx + 1;
        unsigned back_index = idx + numX*numY;
        unsigned back_right_index = idx + 1 + numX*numY;

        std::vector<boost::shared_ptr<Vertex> > poly_vertices;
        poly_vertices.push_back(vertices[front_index]);
        poly_vertices.push_back(vertices[front_right_index]);
        poly_vertices.push_back(vertices[back_right_index]);
        poly_vertices.push_back(vertices[back_index]);
        polygons.push_back(Polygon::Create(poly_vertices));
    }

    // Top face
    for(unsigned idx=0; idx< numX - 1; idx++)
    {
        unsigned front_index = idx + numX*(numY-1);
        unsigned front_right_index = idx + 1 + numX*(numY-1);
        unsigned back_index = idx + + numX*(numY-1) + numX*numY;
        unsigned back_right_index = idx + numX*(numY-1) + 1 + numX*numY;

        std::vector<boost::shared_ptr<Vertex> > poly_vertices;
        poly_vertices.push_back(vertices[front_index]);
        poly_vertices.push_back(vertices[front_right_index]);
        poly_vertices.push_back(vertices[back_right_index]);
        poly_vertices.push_back(vertices[back_index]);
        polygons.push_back(Polygon::Create(poly_vertices));
    }

    // Create a part
    boost::shared_ptr<Part<DIM> > p_part = Part<DIM>::Create();
    for(unsigned idx=0; idx<polygons.size(); idx++)
    {
        p_part->AddPolygon(polygons[idx], true);
    }
    return p_part;
}

template<unsigned DIM>
boost::shared_ptr<Part<DIM> > MappableGridGenerator<DIM>::GenerateCylinder(double cylinder_radius,
                                               double cylinder_thickness,
                                               double cylinder_angle,
                                               double cylinder_height,
                                               unsigned numX,
                                               unsigned numY)
{
    boost::shared_ptr<Part<DIM> > p_part = GeneratePlane(numX, numY);

    // Get the part extents
    c_vector<double, 2*DIM> bbox = p_part->GetBoundingBox();

    // Get the vertices
    std::vector<boost::shared_ptr<Vertex> > vertices = p_part->GetVertices();
    for(unsigned idx =0; idx<vertices.size(); idx++)
    {
        double x_frac = vertices[idx]->rGetLocation()[0] / (bbox[1] - bbox[0]);
        double angle = x_frac * cylinder_angle;

        double y_frac = vertices[idx]->rGetLocation()[1] / (bbox[3] - bbox[2]);
        double height = y_frac * cylinder_height;

        double z_frac = vertices[idx]->rGetLocation()[2] / (bbox[5] - bbox[4]);
        double radius = cylinder_radius - cylinder_thickness * z_frac;

        // Get the new x
        c_vector<double, DIM> new_position;
        new_position[0] = radius * std::cos(angle);

        // Get the new y
        new_position[1] = height;

        // Get the new z
        new_position[2] = radius * std::sin(angle);

        vertices[idx]->Translate(new_position - vertices[idx]->rGetLocation());
    }
    return p_part;
}

template<unsigned DIM>
boost::shared_ptr<Part<DIM> > MappableGridGenerator<DIM>::GenerateHemisphere(double sphere_radius,
                                               double sphere_thickness,
                                               double sphere_azimuth_angle,
                                               double sphere_polar_angle,
                                               unsigned numX,
                                               unsigned numY)
{
    boost::shared_ptr<Part<DIM> > p_part = GeneratePlane(numX, numY);

    // The the part extents
    c_vector<double, 2*DIM> bbox = p_part->GetBoundingBox();

    // Get the vertices
    std::vector<boost::shared_ptr<Vertex> > vertices = p_part->GetVertices();
    for(unsigned idx =0; idx<vertices.size(); idx++)
    {
        double x_frac = vertices[idx]->rGetLocation()[0] / (bbox[1] - bbox[0]);
        double azimuth_angle = x_frac * sphere_azimuth_angle;

        double y_frac = vertices[idx]->rGetLocation()[1] / (bbox[3] - bbox[2]);
        double polar_angle = y_frac * sphere_polar_angle;

        double z_frac = vertices[idx]->rGetLocation()[2] / (bbox[5] - bbox[4]);
        double radius = sphere_radius - sphere_thickness * z_frac;

        // Get the new x
        c_vector<double, DIM> new_position;
        new_position[0] = radius * std::cos(azimuth_angle) * std::sin(polar_angle);

        // Get the new y
        new_position[1] = radius * std::cos(polar_angle);

        // Get the new z
        new_position[2] = radius * std::sin(azimuth_angle) * std::sin(polar_angle);

        vertices[idx]->Translate(new_position - vertices[idx]->rGetLocation());
    }
    return p_part;
}

//Explicit instantiation
template class MappableGridGenerator<2> ;
template class MappableGridGenerator<3> ;