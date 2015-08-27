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
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include "Polygon.hpp"

#include "VoronoiGenerator.hpp"

// Hang Si's tetgen.

#define REAL double
#define VOID void
#include "tetgen15.h"
#undef REAL
#undef VOID

struct triangulateio;

template<unsigned DIM>
VoronoiGenerator<DIM>::VoronoiGenerator()
{
}

template<unsigned DIM>
VoronoiGenerator<DIM>::~VoronoiGenerator()
{

}

template<unsigned DIM>
boost::shared_ptr<Part<DIM> > VoronoiGenerator<DIM>::Generate(boost::shared_ptr<Part<DIM> > pPart,
        std::vector<boost::shared_ptr<Vertex> > seeds, unsigned numSeeds)
{
    if(DIM==2)
    {
        EXCEPTION("This generator works in 3D only");
    }

    boost::shared_ptr<Part<DIM> > p_tesselation = Part<DIM>::Create();

    c_vector<double, 2*DIM> extents = pPart->GetBoundingBox();

    // If no seeds have been provided generate some random ones
    if(seeds.size() == 0)
    {
        boost::mt19937 rng;
        boost::uniform_real<> u_dist(0.0, 1.0);
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > rand_uniform(rng, u_dist);
        for(unsigned idx = 0; idx< numSeeds; idx++)
        {
            c_vector<double, DIM> location;
            for(unsigned jdx=0; jdx<DIM; jdx++)
            {
                location[jdx] = extents[2*jdx] + (extents[2*jdx + 1] - extents[2*jdx])*rand_uniform();
            }
            seeds.push_back(Vertex::Create(location));
        }
    }

    // Use tetgen to make the tessellation
    class tetgen15::tetgenio mesher_input, mesher_output;

    mesher_input.initialize();
    mesher_output.initialize();
    mesher_input.numberofpoints = seeds.size();
    mesher_input.pointlist = new double[mesher_input.numberofpoints * DIM];

    for (unsigned idx=0; idx<seeds.size(); idx++)
    {
        for (unsigned jdx=0; jdx < DIM; jdx++)
        {
            mesher_input.pointlist[DIM*idx + jdx] = seeds[idx]->rGetLocation()[jdx];
        }
    }
    tetgen15::tetrahedralize((char*)"veeQ", &mesher_input, &mesher_output);

    // Create 2-point polygons corresponding to each edge
    std::vector<boost::shared_ptr<Polygon> > polygons;
    for (int n=0; n<mesher_output.numberofvedges; ++n) {
      tetgen15::tetgenio::voroedge e =  mesher_output.vedgelist[n];
      int n0 = e.v1;
      int n1 = e.v2;
      double* u = &mesher_output.vpointlist[DIM*n0];
      double* v = n1 == -1 ? e.vnormal : &mesher_output.vpointlist[DIM*n1]; // -1 indicates ray

      boost::shared_ptr<Vertex> p_vertex1 = Vertex::Create(u[0], u[1], u[2]);
      boost::shared_ptr<Vertex> p_vertex2 = Vertex::Create(v[0], v[1], v[2]);

      // Only include edges that are inside the domain.
      bool vert1_inside = true;
      bool vert2_inside = true;
      if(p_vertex1->rGetLocation()[0] < extents[0] || p_vertex1->rGetLocation()[0] > extents[1])
      {
          vert1_inside = false;
      }
      if(p_vertex1->rGetLocation()[1] < extents[2] || p_vertex1->rGetLocation()[1] > extents[3])
      {
          vert1_inside = false;
      }
      if(p_vertex2->rGetLocation()[0] < extents[0] || p_vertex2->rGetLocation()[0] > extents[1])
      {
          vert2_inside = false;
      }
      if(p_vertex2->rGetLocation()[1] < extents[2] || p_vertex2->rGetLocation()[1] > extents[3])
      {
          vert2_inside = false;
      }
      if(DIM==3)
      {
          if(p_vertex1->rGetLocation()[2] < extents[4] || p_vertex1->rGetLocation()[2] > extents[5])
          {
              vert1_inside = false;
          }
          if(p_vertex2->rGetLocation()[2] < extents[4] || p_vertex2->rGetLocation()[2] > extents[5])
          {
              vert2_inside = false;
          }
      }

      if(vert1_inside && vert2_inside)
      {
          if(!vert1_inside)
          {
              if(p_vertex1->rGetLocation()[0] < extents[0])
              {
                  p_vertex1->SetCoordinate(0, extents[0]);
              }
              if(p_vertex1->rGetLocation()[0] > extents[1])
              {
                  p_vertex1->SetCoordinate(0, extents[1]);
              }
              if(p_vertex1->rGetLocation()[1] < extents[2])
              {
                  p_vertex1->SetCoordinate(1, extents[2]);
              }
              if(p_vertex1->rGetLocation()[1] > extents[3])
              {
                  p_vertex1->SetCoordinate(1, extents[3]);
              }
              if(DIM==3)
              {
                  if(p_vertex1->rGetLocation()[2] < extents[4])
                  {
                      p_vertex1->SetCoordinate(2, extents[4]);
                  }
                  if(p_vertex1->rGetLocation()[2] > extents[5])
                  {
                      p_vertex1->SetCoordinate(2, extents[5]);
                  }
              }

          }
          if(!vert2_inside)
          {
              if(p_vertex2->rGetLocation()[0] < extents[0])
              {
                  p_vertex2->SetCoordinate(0, extents[0]);
              }
              if(p_vertex2->rGetLocation()[0] > extents[1])
              {
                  p_vertex2->SetCoordinate(0, extents[1]);
              }
              if(p_vertex2->rGetLocation()[1] < extents[2])
              {
                  p_vertex2->SetCoordinate(1, extents[2]);
              }
              if(p_vertex2->rGetLocation()[1] > extents[3])
              {
                  p_vertex2->SetCoordinate(1, extents[3]);
              }
              if(DIM==3)
              {
                  if(p_vertex2->rGetLocation()[2] < extents[4])
                  {
                      p_vertex2->SetCoordinate(2, extents[4]);
                  }
                  if(p_vertex2->rGetLocation()[2] > extents[5])
                  {
                      p_vertex2->SetCoordinate(2, extents[5]);
                  }
              }

          }

          boost::shared_ptr<Polygon> p_polygon = Polygon::Create(p_vertex1);
          p_polygon->AddVertex(p_vertex2);
          polygons.push_back(p_polygon);
      }
    }

    // Add the polygons to the part
    for (unsigned idx=0; idx<polygons.size(); idx++)
    {
        p_tesselation->AddPolygon(polygons[idx]);
    }
    return p_tesselation;
}

//Explicit instantiation
template class VoronoiGenerator<2> ;
template class VoronoiGenerator<3> ;