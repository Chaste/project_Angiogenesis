#!/usr/bin/env python

"""Copyright (c) 2005-2016, University of Oxford.
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

.. module:: image_to_surface
    :synopsis: Tools for extracting geometric features from images/pixel based features
"""
import sys
import vtk
from vmtk import vmtkscripts
import chaste.utility.bases as bases
import chaste.interfaces.vtk_tools

class VtkImageToPolyData2d(bases.SimpleIOBase):
    """Convert a 2D image in vtk format to a boundary in polydata format.
    
    The boundary delineates the regions with pixel values above and below the specified threshold. There are
    also options for controlling the quality of the extracted boundary by successive smoothing, subdivision 
    and decimation. Suitable choices of these parameters may require some trial and error.
    """
    def __init__(self):
        """Convert a 2D image in vtk format to a boundary in polydata format.
        
        The image is thresholded, subject to decimation and then linear subdivision. Changing paramters for
        the latter two filters can change the quality of the final boundary description.
        
        self.input is 2D vtk image data
        self.output is a collection of polylines as vtk polydata containing boundaries
        self.target_reduction is the fraction of edges to remove when decimating the boundary
        self.num_subdivisions is the number of linear subdivisions
        self.resize_factors [xyz] factors to scale the image by before processing
        self.clipping_factor fraction of the original image to remove, useful to avoid edge effects
        self.threshold the threshold value across which the boundary lies
        self.smoothing_length can be tuned to smooth the boundary after extraction, too high and edges are rounded
        """
        super(VtkImageToPolyData2d, self).__init__()
        self.tool_name = "Chaste"+self.__class__.__name__
        
        self.target_reduction = 0.9994
        self.num_subdivisions = 3
        self.resize_factors = None
        self.clipping_factor = 0.0
        self.threshold = 1.0
        self.smoothing_length = 3.0
    
    def update(self):
        """Run the tool"""
        # Size the image if needed
        if self.resize_factors is not None:
            reduction_filter = vtk.vtkImageResize()
            reduction_filter.SetInput(self.input)
            reduction_filter.SetMagnificationFactors(self.resize_factors[0], self.resize_factors[1], self.resize_factors[2])
            reduction_filter.SetResizeMethodToMagnificationFactors()
            reduction_filter.Update()
        
        # Threshold
        threshold = vtk.vtkThreshold()
        if self.resize_factors is not None:
            threshold.SetInputConnection(reduction_filter.GetOutputPort())
        else:
            threshold.SetInput(self.input)
        threshold.ThresholdByLower(1.0)
        threshold.Update()
        
        # Convert to polydata
        surface = vtk.vtkGeometryFilter()
        surface.SetInputConnection(threshold.GetOutputPort())
        surface.Update()
        
        # Triangulate
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInputConnection(surface.GetOutputPort())
        triangle.Update()
    
        # Decimate
        decimate = vtk.vtkDecimatePro()
        decimate.SetInputConnection(triangle.GetOutputPort())
        decimate.SetTargetReduction(self.target_reduction)
        decimate.SetFeatureAngle(15.0)
        decimate.Update()
        
        # Do loop subdivision
        su = vtk.vtkLinearSubdivisionFilter()
        su.SetInputConnection(decimate.GetOutputPort())
        su.SetNumberOfSubdivisions(self.num_subdivisions)
        su.Update()
        
        # Clip the boundaries, recommended to ensure straight inlets and outlets
        if self.clipping_factor > 0.0:
            bounds = su.GetOutput().GetBounds()
            width = bounds[1] - bounds[0]
            height = bounds[3] - bounds[2]
            
            p1 = vtk.vtkPlane()
            p1.SetOrigin(bounds[0] + self.clipping_factor*width, 0.0, 0)
            p1.SetNormal(1, 0, 0)
            p2 = vtk.vtkPlane()
            p2.SetOrigin(0.0, bounds[2] + self.clipping_factor*height, 0)
            p2.SetNormal(0, 1, 0)            
            p3 = vtk.vtkPlane()
            p3.SetOrigin(bounds[1] - self.clipping_factor*width, 0.0, 0)
            p3.SetNormal(-1, 0, 0)        
            p4 = vtk.vtkPlane()
            p4.SetOrigin(0.0, bounds[3] - self.clipping_factor*height, 0)
            p4.SetNormal(0, -1, 0)        
        
            c1 = vtk.vtkClipPolyData()
            c1.SetInputConnection(su.GetOutputPort())
            c1.SetClipFunction(p1)
            c1.SetValue(0.0)            
            c2 = vtk.vtkClipPolyData()
            c2.SetInputConnection(c1.GetOutputPort())
            c2.SetClipFunction(p2)
            c2.SetValue(0.0) 
            c3 = vtk.vtkClipPolyData()
            c3.SetInputConnection(c2.GetOutputPort())
            c3.SetClipFunction(p3)
            c3.SetValue(0.0)
            c4 = vtk.vtkClipPolyData()
            c4.SetInputConnection(c3.GetOutputPort())
            c4.SetClipFunction(p4)
            c4.SetValue(0.0)
            
            su = vtk.vtkGeometryFilter()
            su.SetInputConnection(c4.GetOutputPort())
            su.Update()
        
        # Get the boundary as triangulated polydata
        feature_edges = vtk.vtkFeatureEdges()
        feature_edges.SetInputConnection(su.GetOutputPort())  
        feature_edges.Update()   
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(feature_edges.GetOutputPort())
        clean.Update()
         
        triangle2 = vtk.vtkTriangleFilter()
        triangle2.SetInputConnection(clean.GetOutputPort())
        triangle2.Update()
        self.output = triangle2.GetOutput()
        
        # Convert the boundary to a polyline description and smooth it
        self.output = chaste.interfaces.vtk_tools.vtk_tools.vtk_lines_to_polylines_region(self.output)
        smoother = vmtkscripts.vmtkCenterlineResampling()
        smoother.Centerlines = self.output
        smoother.length = self.smoothing_length
        smoother.Execute()
        self.output = smoother.Centerlines
        
if __name__=='__main__':
    
    tool = VtkImageToPolyData2d()
    tool.run_standalone(sys.argv)
    