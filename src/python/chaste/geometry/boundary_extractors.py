import vtk
    
class BoundaryExtractor2d():
    
    
    def __init__(self):
        
        self.labels = None
        
    def update(self):
        
        """
        For each marked location on the boudaries extract the start and end point of the 
        boundary edge.
        """
        
        self.output = []
        
        for eachLabel in self.labels:
            threshold = vtk.vtkThreshold()
            threshold.SetInput(self.input)
            threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_POINTS", "PointBoundaryLabel")
            threshold.ThresholdBetween(eachLabel, eachLabel)
            threshold.Update()
            
            connectivity = vtk.vtkConnectivityFilter()
            connectivity.SetInputConnection(threshold.GetOutputPort())
            connectivity.SetExtractionModeToAllRegions()
            connectivity.ColorRegionsOn()
            connectivity.Update()
             
            polysurface = vtk.vtkGeometryFilter()
            polysurface.SetInputConnection(connectivity.GetOutputPort())
            polysurface.Update()
             
            points = polysurface.GetOutput().GetPoints()
            num_points = polysurface.GetOutput().GetNumberOfPoints()   
            point_label_array = polysurface.GetOutput().GetPointData().GetArray("RegionId")
            point_locs = []
            point_connectivity = []
            point_labels = []
             
            for i in range(num_points):
                point_locs.append([points.GetPoint(i)[0], points.GetPoint(i)[1]])
                point_connectivity.append(0)
                point_labels.append(int(point_label_array.GetTuple1(i)))
                        
            numCells = polysurface.GetOutput().GetNumberOfLines()  
            cellArray = polysurface.GetOutput().GetLines()
            cellArray.InitTraversal()
            segList = vtk.vtkIdList()
                 
            for i in range(numCells): 
                cellArray.GetNextCell(segList)
                for j in range(0, segList.GetNumberOfIds()):
                    seg_id = segList.GetId(j)
                    point_connectivity[int(seg_id)] += 1
                     
            # Get point ids
            point_set = list(set(point_labels))
            region_points = []
            for eachRegion in point_set:
                my_points = []
                for idx, eachPoint in enumerate(point_labels):
                    if eachPoint == eachRegion:
                        if point_connectivity[idx] == 1:
                            my_points.append(point_locs[idx])
                region_points.append(my_points)
                 
            self.output.append(region_points)
    
class BoundaryExtract():
    
    def __init__(self):

        pass
    
    def generate_lines(self, file_path, label):
        
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file_path)
        reader.Update()
        
        threshold = vtk.vtkThreshold()
        threshold.SetInputConnection(reader.GetOutputPort())
        threshold.ThresholdBetween(label, label)
        threshold.Update()
        
        connectivity = vtk.vtkConnectivityFilter()
        connectivity.SetInputConnection(threshold.GetOutputPort())
        connectivity.SetExtractionModeToAllRegions()
        connectivity.ColorRegionsOn()
        connectivity.Update()
        
        surface = vtk.vtkGeometryFilter()
        surface.SetInputConnection(connectivity.GetOutputPort())
        surface.Update()
        
        points = surface.GetOutput().GetPoints()
        num_points = surface.GetOutput().GetNumberOfPoints()   
        point_label_array = surface.GetOutput().GetPointData().GetArray("RegionId")
        point_locs = []
        point_connectivity = []
        point_labels = []
        
        for i in range(num_points):
            point_locs.append([points.GetPoint(i)[0], points.GetPoint(i)[1]])
            point_connectivity.append(0)
            point_labels.append(int(point_label_array.GetTuple1(i)))
                   
        numCells = surface.GetOutput().GetNumberOfLines()  
        cellArray = surface.GetOutput().GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
            
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_connectivity[int(seg_id)] += 1
                
        # Get point ids
        point_set = list(set(point_labels))
        region_points = []
        for eachRegion in point_set:
            my_points = []
            for idx, eachPoint in enumerate(point_labels):
                if eachPoint == eachRegion:
                    if point_connectivity[idx] == 1:
                        my_points.append(point_locs[idx])
            region_points.append(my_points)
            
        return region_points
    
class BoundaryExtract3d():
    
    def __init__(self):

        self.surface = None
        self.regions = None
        self.labels = None
    
    def update(self):
        
        self.regions = []
        for eachLabel in self.labels :
            threshold = vtk.vtkThreshold()
            threshold.SetInput(self.surface)
            threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_CELLS", "CellEntityIds")
            threshold.ThresholdBetween(eachLabel, eachLabel)
            threshold.Update()
        
            surface = vtk.vtkGeometryFilter()
            surface.SetInputConnection(threshold.GetOutputPort())
            surface.Update()
            
            # Triangulate
            triangle = vtk.vtkTriangleFilter()
            triangle.SetInputConnection(surface.GetOutputPort())
            triangle.Update()
        
            su = vtk.vtkLinearSubdivisionFilter()
            su.SetInputConnection(triangle.GetOutputPort())
            su.SetNumberOfSubdivisions(3)
            su.Update()
            
            self.regions.append(su.GetOutput())
            
    def get_output(self):
        
        return self.regions