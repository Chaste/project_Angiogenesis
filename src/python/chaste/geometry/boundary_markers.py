import numpy as np
import vtk
import chaste.mesh.converters
import chaste.utility.bases as bases

class BoundaryMarker2d(bases.SimpleIOBase):
    
    def __init__(self):
        super(BoundaryMarker2d, self).__init__()
        self.boundary_edges = None
        self.seed_points = None
        self.feature_angle = (10.0 /180.0)*np.pi
        self.labels = range(1,11)
    
    def unit_vector(self, vector):
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def crawl(self, input_point, previous_point, connectivity, edges):
        cell_1 = connectivity[input_point][0]
        cell_2 = connectivity[input_point][1]
        if edges[cell_1][0] == input_point and edges[cell_1][1] != previous_point:
            opposite = edges[cell_1][1]
            opposite_cell = cell_1
        elif edges[cell_1][1] == input_point and edges[cell_1][0] != previous_point:
            opposite = edges[cell_1][0]
            opposite_cell = cell_1
        elif edges[cell_2][0] == input_point and edges[cell_2][1] != previous_point:
            opposite = edges[cell_2][1]
            opposite_cell = cell_2        
        elif edges[cell_2][1] == input_point and edges[cell_2][0] != previous_point:
            opposite = edges[cell_2][0]
            opposite_cell = cell_2                 
        return opposite, opposite_cell

    def update(self):
        
        """
        For each location in the point_seeds input find the closest pair of edges on the surface, then
        mark all edges at an angle less than the feature angle to the chosen edges. Marker labels
        can be optionally specified.
        """
        
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInput(self.input)
        triangle.Update()
        line_boundary = triangle.GetOutput()
        
        converter = chaste.mesh.converters.VtkToTriMesh()
        converter.input = line_boundary
        points, edges = converter.update()
        
        # build point edge connectivity
        connectivity = []
        for idx in range(len(points)):
            connectivity.append([])
        
        cell_markers = []
        for idx, eachEdge in enumerate(edges):
            cell_markers.append(0)
            connectivity[eachEdge[0]].append(idx)
            connectivity[eachEdge[1]].append(idx)
        
        seed_neighbour_ids = []
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(self.input)
        for eachLoc in self.seed_points:
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            seed_neighbour_ids.append(closest_id)
        
        # Mark planar cells for each input point
        self.boundary_edges = []
        interior_points = np.zeros(len(points))
        boundary_points = np.zeros(len(points))
        
        for idx, eachStartPoint in enumerate(seed_neighbour_ids):
            
            if idx<len(self.labels):
                this_label = self.labels[idx]
            else:
                this_label = len(self.labels)
            
            cell_1 = connectivity[eachStartPoint][0]
            cell_2 = connectivity[eachStartPoint][1]
            
            if edges[cell_1][0] == eachStartPoint:
                opposite1 = edges[cell_1][1]
            else:
                opposite1 = edges[cell_1][0]
                
            if edges[cell_2][0] == eachStartPoint:
                opposite2 = edges[cell_2][1]
            else:
                opposite2 = edges[cell_2][0]
            dir1 = self.unit_vector(np.array(points[eachStartPoint]) - np.array(points[opposite1]))
            dir2 = self.unit_vector(np.array(points[eachStartPoint]) - np.array(points[opposite2]))
            input_angle = self.angle_between(dir1, dir2)
            
            if abs(input_angle-np.pi)<=self.feature_angle:
                cell_markers[cell_1] = this_label
                cell_markers[cell_2] = this_label
                
                interior_points[eachStartPoint] = 1
                boundary_points[eachStartPoint] = this_label
                input_point = opposite1
                previous_point = eachStartPoint
                previous_dir = dir1
                edge_found = False
                
                while not edge_found:
                    opposite, opposite_cell = self.crawl(input_point, previous_point, connectivity, edges)
                    opposite_dir = self.unit_vector(np.array(points[input_point]) - np.array(points[opposite]))
                    new_angle = self.angle_between(previous_dir, opposite_dir)
                    
                    if abs(new_angle-np.pi)<=self.feature_angle or abs(new_angle)<=self.feature_angle:
                        cell_markers[opposite_cell] = this_label
                        previous_point = input_point
                        input_point = opposite
                        previous_dir = opposite_dir
                        interior_points[previous_point] = 1
                        boundary_points[previous_point] = this_label
                    else:
                        boundary_points[input_point] = this_label
                        edge_found = True
                        right_point = input_point 
                        break
                    
                input_point = opposite2
                previous_point = eachStartPoint
                previous_dir = dir2
                edge_found = False
                
                while not edge_found:
                    opposite, opposite_cell = self.crawl(input_point, previous_point, connectivity, edges)
                    opposite_dir = self.unit_vector(np.array(points[input_point]) - np.array(points[opposite]))
                    new_angle = self.angle_between(previous_dir, opposite_dir)
                    
                    if abs(new_angle-np.pi)<=self.feature_angle or abs(new_angle)<=self.feature_angle:
                        cell_markers[opposite_cell] = this_label
                        previous_point = input_point
                        input_point = opposite
                        previous_dir = opposite_dir
                        interior_points[previous_point] = 1
                        boundary_points[previous_point] = this_label
                    else:
                        edge_found = True
                        left_point = input_point
                        boundary_points[input_point] = this_label
                        break  
                    
                self.boundary_edges.append([(left_point, right_point), this_label]) 
            else:
                interior_points[eachStartPoint] = 1
                boundary_points[eachStartPoint] = this_label
                
                # Get the left and right cell centres
                mid1 = (np.array(points[eachStartPoint]) + np.array(points[opposite1]))/2.0
                mid2 = (np.array(points[eachStartPoint]) + np.array(points[opposite2]))/2.0
                dist1 = np.linalg.norm(np.array(self.seed_points[idx]) - mid1)
                dist2 = np.linalg.norm(np.array(self.seed_points[idx]) - mid2)
                if dist1 < dist2:
                    cell_markers[cell_1] = this_label
                    left_point = opposite1
                    right_point = eachStartPoint
                    boundary_points[opposite1] = this_label
                else:
                    cell_markers[cell_2] = this_label
                    left_point = opposite2
                    right_point = eachStartPoint
                    boundary_points[opposite2] = this_label
                left_point = opposite2
                right_point = opposite1 
                self.boundary_edges.append([(left_point, right_point), this_label]) 
            
        # Mark the vtk_data using cell markers
        boundary_label = vtk.vtkFloatArray()
        boundary_label.SetNumberOfComponents(1)
        boundary_label.SetName("CellBoundaryLabel")
        
        interior_label = vtk.vtkFloatArray()
        interior_label.SetNumberOfComponents(1)
        interior_label.SetName("InteriorLabel")
        
        boundary_point_label = vtk.vtkFloatArray()
        boundary_point_label.SetNumberOfComponents(1)
        boundary_point_label.SetName("PointBoundaryLabel")
        
        for idx, eachEdge in enumerate(edges):
            boundary_label.InsertNextTupleValue((float(cell_markers[idx]),))  
            
        for idx in range(len(interior_points)):
            interior_label.InsertNextTupleValue((float(interior_points[idx]),))  
            
        for idx in range(len(boundary_points)):
            boundary_point_label.InsertNextTupleValue((float(boundary_points[idx]),))  
        
        
        line_boundary.GetCellData().AddArray(boundary_label) 
        line_boundary.GetPointData().SetScalars(interior_label) 
        line_boundary.GetPointData().AddArray(boundary_point_label) 
        self.output = line_boundary

class BoundaryMarker3d():
    
    def __init__(self, mesh, domains=None, domain_labels=None, boundarys=None, boundary_labels = None):
        self.mesh = mesh
        self.domains = domains
        self.domain_labels = domain_labels
        self.boundarys = boundarys
        self.boundary_labels = boundary_labels
        
    def generate(self):
        
        # Mark the domains
        mesh_func = df.MeshFunction("size_t", self.mesh, 3)
        
#         if len(self.element_labels)>0:
#             num_cells = len(self.mesh.cells())
#             for idx in range(num_cells):
#                 mesh_func[idx] = int(self.element_labels[idx])
    
        # Set up facet markers
        boundaries = df.FacetFunction("size_t", self.mesh)
        boundaries.set_all(0)
        
        if self.boundarys is not None:
            for idx, eachBoundary in enumerate(self.boundarys):
                eachBoundary.mark(boundaries, self.boundary_labels[idx])
                
        edges = df.EdgeFunction("size_t", self.mesh)
        edges.set_all(0)

        return mesh_func, boundaries, edges