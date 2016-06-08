import vtk
import numpy as np
from vmtk import vmtkscripts, pypes
import chaste.utility.bases as bases
import chaste.utility.readwrite

class VtkImageToPolyData2d(bases.SimpleIOBase):
    
    """ Convert an image in vtk format to polydata by extracting the boundary between regions on
    either side of a specified threshold. Optional steps include boundary smoothing, resizing and ROI extraction.
    """
    
    def __init__(self):
        
        """
        Input is 2D Vtk image data
        Output is vtk polydata with the boundary between the regions separated by the specified threshold
        """
        
        super(VtkImageToPolyData2d, self).__init__()
        self.surface = None
        self.target_reduction = 0.9994
        self.num_subdivisions = 3
        self.resize_factors = None
        self.clipping_factor = 0.0
        self.threshold = 1.0
    
    def update(self):
        
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
# 
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
        
        feature_edges = vtk.vtkFeatureEdges()
        feature_edges.SetInputConnection(su.GetOutputPort())  
        feature_edges.Update()   
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(feature_edges.GetOutputPort())
        clean.Update()
         
        triangle2 = vtk.vtkTriangleFilter()
        triangle2.SetInputConnection(clean.GetOutputPort())
        triangle2.Update()
        
        self.surface = su.GetOutput()
        self.output = triangle2.GetOutput()
        self.output = self.boundary_to_polylines()
        
        smoother = vmtkscripts.vmtkCenterlineResampling()
        smoother.Centerlines = self.output
        smoother.length = 3.0
        smoother.Execute()
        self.output = smoother.Centerlines
        
    def boundary_to_polylines(self):
        
        """
        Convert the shape boundary from a collection of vtk lines to 
        a collection of vtk polylines. This allows for subsequent smoothing of the boundary.
        """
        
        numPoints = self.output.GetNumberOfPoints()    
        vtkpoints = self.output.GetPoints()  
        connectivity = []
        
        points = [] 
        for i in range(numPoints):
            points.append(vtkpoints.GetPoint(i))
            connectivity.append([])
            
        numCells = self.output.GetNumberOfLines()  
        cellArray = self.output.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
            
        edges = []
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(int(seg_id))
            edges.append((point_indices[0], point_indices[1]))
            connectivity[point_indices[0]].append(i)
            connectivity[point_indices[1]].append(i)
            
        regions = []
        point_visited = np.zeros(len(points))
        
        for idx, eachPoint in enumerate(points):
            if point_visited[idx] == 0:
                point_visited[idx] = 1
                region_points = [idx]
                current_point_id = idx
                previous_point_id = idx
                found_visited = False
                while not found_visited:
                    edge_1 = edges[connectivity[current_point_id][0]]
                    edge_2 = edges[connectivity[current_point_id][1]]
                    if edge_1[0] == current_point_id:
                        opp1 = edge_1[1]
                    else:
                        opp1 = edge_1[0]
                    if opp1 != previous_point_id:
                        next_edge = edge_1
                    else:
                        next_edge = edge_2
                              
                    if next_edge[0]== current_point_id:
                        next_point_id = next_edge[1]  
                    else:
                        next_point_id = next_edge[0]     
                              
                    if point_visited[next_point_id] == 1:
                        region_points.append(next_point_id)
                        break
                      
                    point_visited[next_point_id] = 1
                    region_points.append(next_point_id)
                    previous_point_id = current_point_id
                    current_point_id = next_point_id
                          
                if len(region_points) > 1:
                    regions.append(region_points)
                        
        new_points = vtk.vtkPoints()
        new_points.SetNumberOfPoints(len(points))
        for idx, eachPoint in enumerate(points):
            new_points.SetPoint(idx, points[idx][0], points[idx][1], 0.0)
            
        lines = vtk.vtkCellArray()
        for eachRegion in regions:
            lines.InsertNextCell(len(eachRegion))
            for eachPoint in eachRegion:
                lines.InsertCellPoint(eachPoint)
                
        polygon = vtk.vtkPolyData()
        polygon.SetPoints(new_points)
        polygon.SetLines(lines)
        
        return polygon
    
class VtkImageToPolyData3d(bases.SimpleIOBase):
    
    def __init__(self):
        super(VtkImageToPolyData2d, self).__init__()
        self.temp_dir = "."
        self.pad_size = 5
    
    def update(self):     
        
        if self.pad_size > 0:
            pad_filter = vtk.vtkImageConstantPad()
            pad_filter.SetInput(self.input)
            pad_filter.SetConstant(0)
            pad_filter.SetOutputWholeExtent(self.input.GetExtent()[0]-self.pad_size, 
                                            self.input.GetExtent()[1]+self.pad_size, 
                                            self.input.GetExtent()[2]-self.pad_size, 
                                            self.input.GetExtent()[3]+self.pad_size, 
                                            self.input.GetExtent()[4]-self.pad_size, 
                                            self.input.GetExtent()[5]+self.pad_size)
            pad_filter.Update()
            self.input = pad_filter.GetOutput()
        
        chaste.utility.readwrite.write(self.input, self.temp_dir+"/temp_input_image.vti")
         
        myArguments = 'vmtklevelsetsegmentation -ifile ' + self.temp_dir + "/temp_input_image.vti" + ' -ofile ' + self.temp_dir +"/temp_level_set.vti"
        pypes.PypeRun(myArguments)
        myArguments = 'vmtkmarchingcubes -ifile ' + self.temp_dir +"/temp_level_set.vti" + ' -ofile ' + self.temp_dir + "/temp_surface.vtp"
        pypes.PypeRun(myArguments)
        
        self.output = chaste.utility.readwrite.read(self.temp_dir + "/temp_surface.vtp")
        return self.output