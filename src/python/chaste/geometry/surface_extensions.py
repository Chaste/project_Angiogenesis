import numpy as np
import vtk

import chaste.geometry.converters.other
import chaste.geometry.boundary_markers

class FlowExtension2d():
    
    def __init__(self):
        
        self.length_factor = 1.5
        self.division_length = 4.0
        self.open_surface = None
        self.extended_surface = None
        self.reference_surface = None
        self.edges = None
        self.midpoints = None
    
    def update(self):
        
        # create a new vtk object
        self.extended_surface = vtk.vtkPolyData() 
        
        converter = chaste.geometry.other.VtkToTri()
        old_points, old_edges = converter.generate(self.open_surface)
        
        bounds = self.open_surface.GetBounds()
        mid_x = (bounds[1] - bounds[0])/2.0 + bounds[0]
        mid_y = (bounds[3] - bounds[2])/2.0 + bounds[2]
        
        new_vtk_points = vtk.vtkPoints()
        new_vtk_lines = vtk.vtkCellArray()
        
        for eachEdge in old_edges:
            line = vtk.vtkLine()
            loc1 = (old_points[eachEdge[0]][0], old_points[eachEdge[0]][1], 0.0)
            loc2 = (old_points[eachEdge[1]][0], old_points[eachEdge[1]][1], 0.0)
            pointId1 = new_vtk_points.InsertNextPoint(loc1)
            pointId2 = new_vtk_points.InsertNextPoint(loc2)
            line.GetPointIds().InsertId(0, pointId1)
            line.GetPointIds().InsertId(1, pointId2)
            new_vtk_lines.InsertNextCell(line)
                
        self.extended_surface.SetPoints(new_vtk_points)
        self.extended_surface.SetLines(new_vtk_lines)
        
        locator = vtk.vtkCellLocator()
        locator.SetDataSet(self.reference_surface)
        locator.BuildLocator()
 
        self.midpoints = []
        self.new_inlet_points= []
        self.new_outlet_points= []
        
        for region_index, eachBoundaryLabel in enumerate(self.edges):
            for eachBoundary in eachBoundaryLabel:
                edge_start_point = np.array(eachBoundary[0])
                edge_end_point = np.array(eachBoundary[1])
                edge_length = np.linalg.norm(edge_start_point - edge_end_point)
                edge_tangent = (edge_start_point - edge_end_point) / edge_length
                edge_normal = np.array((edge_tangent[1], -edge_tangent[0]))# * self.norm_list[mydex]
                
                vec_mid = np.array((mid_x, mid_y))
                edge_mid = (edge_start_point + edge_end_point)/2.0
                test_point1 = 0.1*edge_length*edge_normal +  edge_mid
                if np.linalg.norm(test_point1 - vec_mid) < np.linalg.norm(edge_mid - vec_mid):
                    edge_normal *= -1.0
                
                # Draw a box along the edge normal
                num_points = int(edge_length * self.length_factor/self.division_length) + 1
                line_1_points = []
                line_2_points = []
            
                for idx in range(num_points):
                    p1 = np.copy(edge_start_point) + self.division_length * edge_normal * float(idx + 1)
                    p2 = np.copy(edge_end_point) + self.division_length * edge_normal * float(idx + 1)
                    line_1_points.append(p1)
                    line_2_points.append(p2)
                    if idx == num_points -1:
                        if region_index + 1 == 1:
                            self.new_inlet_points.append((p1 + p2)/2.0)
                        if region_index + 1 == 2:
                            self.new_outlet_points.append((p1 + p2)/2.0)
                        self.midpoints.append((p1 + p2)/2.0)
                
                # Get the vtk ids of the start and end points
                locator = vtk.vtkKdTreePointLocator()
                locator.SetDataSet(self.extended_surface)
                probe_loc = np.array((edge_start_point[0], edge_start_point[1], 0.0))
                pt_id_1 = locator.FindClosestPoint(probe_loc)
                probe_loc = np.array((edge_end_point[0], edge_end_point[1], 0.0))
                pt_id_2 = locator.FindClosestPoint(probe_loc)
                
                vtk_numPoints = self.extended_surface.GetNumberOfPoints()    
                vtk_points = self.extended_surface.GetPoints()  
                cellArray = self.extended_surface.GetLines()
                
                for idx, eachPoint in enumerate(line_1_points):
                    loc = np.array((eachPoint[0], eachPoint[1], 0.0))
                    vtk_points.InsertNextPoint(loc)
                     
                    line = vtk.vtkLine()  
                    if idx == 0:
                        line.GetPointIds().InsertId(0, pt_id_1)
                    else:
                        line.GetPointIds().InsertId(0, idx + vtk_numPoints-1)
                    line.GetPointIds().InsertId(1, idx + vtk_numPoints)
                    cellArray.InsertNextCell(line)
                     
                for idx, eachPoint in enumerate(line_2_points):
                    loc = np.array((eachPoint[0], eachPoint[1], 0.0))
                    vtk_points.InsertNextPoint(loc)
                     
                    line = vtk.vtkLine()  
                    if idx == 0:
                        line.GetPointIds().InsertId(0, pt_id_2)
                    else:
                        line.GetPointIds().InsertId(0, idx + vtk_numPoints-1 + len(line_1_points))
                    line.GetPointIds().InsertId(1, idx + vtk_numPoints + len(line_1_points))
                    cellArray.InsertNextCell(line)            
                 
                line = vtk.vtkLine()  
                line.GetPointIds().InsertId(0, vtk_numPoints + len(line_1_points) - 1)
                line.GetPointIds().InsertId(1, vtk_numPoints + len(line_1_points) + len(line_2_points) - 1)
                cellArray.InsertNextCell(line) 
                
        # re_label the boundaries
        marker = chaste.geometry.boundary_markers.BoundaryMarker2d()
        marker.surface = self.extended_surface
        marker.inlet_points = self.new_inlet_points
        marker.outlet_points = self.new_outlet_points
        marker.update()
        self.extended_surface = marker.get_output()
                 
    def get_output(self):
        return self.extended_surface
        