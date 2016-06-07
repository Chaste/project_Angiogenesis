import numpy as np
import vtk
import chaste.geometry.converters.other

class BoundaryMarker2d():
    
    def __init__(self):
        self.surface = None
        self.boundary_edges = None
        self.inlet_points = None
        self.outlet_points = None
        self.feature_angle = (10.0 /180.0)*np.pi
    
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
    
    def set_surface(self, surface):
        
        self.surface = surface

    def update(self):
        
        """
        For each location in the point_seeds input find the closest pair of edges on the surface, then
        mark all edges at an angle less than the feature angle to the chosen edges. Marker labels
        can be optionally specified.
        """
        
        converter = chaste.geometry.other.VtkToTri()
        points, edges = converter.generate(self.surface)
        print len(points), len(edges)
        
        # build point edge connectivity
        connectivity = []
        for idx in range(len(points)):
            connectivity.append([])
        
        cell_markers = []
        for idx, eachEdge in enumerate(edges):
            cell_markers.append(0)
            connectivity[eachEdge[0]].append(idx)
            connectivity[eachEdge[1]].append(idx)
            
        print len(connectivity)
        
        inlet_pt_ids = []
        total_ids = []
        original_pts = []
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(self.surface)
        for eachLoc in self.inlet_points:
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            inlet_pt_ids.append(closest_id)
            total_ids.append(closest_id)
            original_pts.append(eachLoc)
        
        outlet_pt_ids = []
        for eachLoc in self.outlet_points:
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            outlet_pt_ids.append(closest_id)
            total_ids.append(closest_id)
            original_pts.append(eachLoc)
        
        # Mark planar cells for each input point
        self.boundary_edges = []
        interior_points = np.zeros(len(points))
        boundary_points = np.zeros(len(points))
        
        for idx, eachStartPoint in enumerate(total_ids):
            
            if eachStartPoint in inlet_pt_ids:
                bound_label = 1
            else:
                bound_label = 2
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
                cell_markers[cell_1] = bound_label
                cell_markers[cell_2] = bound_label
                
                interior_points[eachStartPoint] = 1
                boundary_points[eachStartPoint] = bound_label
                input_point = opposite1
                previous_point = eachStartPoint
                previous_dir = dir1
                edge_found = False
                
                while not edge_found:
                    opposite, opposite_cell = self.crawl(input_point, previous_point, connectivity, edges)
                    opposite_dir = self.unit_vector(np.array(points[input_point]) - np.array(points[opposite]))
                    new_angle = self.angle_between(previous_dir, opposite_dir)
                    
                    if abs(new_angle-np.pi)<=self.feature_angle or abs(new_angle)<=self.feature_angle:
                        cell_markers[opposite_cell] = bound_label
                        previous_point = input_point
                        input_point = opposite
                        previous_dir = opposite_dir
                        interior_points[previous_point] = 1
                        boundary_points[previous_point] = bound_label
                    else:
                        boundary_points[input_point] = bound_label
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
                        cell_markers[opposite_cell] = bound_label
                        previous_point = input_point
                        input_point = opposite
                        previous_dir = opposite_dir
                        interior_points[previous_point] = 1
                        boundary_points[previous_point] = bound_label
                    else:
                        edge_found = True
                        left_point = input_point
                        boundary_points[input_point] = bound_label
                        break  
                    
                self.boundary_edges.append([(left_point, right_point), bound_label]) 
            else:
                interior_points[eachStartPoint] = 1
                boundary_points[eachStartPoint] = bound_label
                
                # Get the left and right cell centres
                mid1 = (np.array(points[eachStartPoint]) + np.array(points[opposite1]))/2.0
                mid2 = (np.array(points[eachStartPoint]) + np.array(points[opposite2]))/2.0
                dist1 = np.linalg.norm(np.array(original_pts[idx]) - mid1)
                dist2 = np.linalg.norm(np.array(original_pts[idx]) - mid2)
                if dist1 < dist2:
                    cell_markers[cell_1] = bound_label
                    left_point = opposite1
                    right_point = eachStartPoint
                    boundary_points[opposite1] = bound_label
                else:
                    cell_markers[cell_2] = bound_label
                    left_point = opposite2
                    right_point = eachStartPoint
                    boundary_points[opposite2] = bound_label
                
                left_point = opposite2
                right_point = opposite1 
                self.boundary_edges.append([(left_point, right_point), bound_label]) 
            
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
            
        for idx in range(interior_points):
            interior_label.InsertNextTupleValue((float(interior_points[idx]),))  
            
        for idx in range(boundary_points):
            boundary_point_label.InsertNextTupleValue((float(boundary_points[idx]),))  
        
        self.surface.GetCellData().AddArray(boundary_label) 
        self.surface.GetPointData().SetScalars(interior_label) 
        self.surface.GetPointData().AddArray(boundary_point_label) 
        self.get_open_surface()
        
    def get_open_surface(self):
        
        threshold = vtk.vtkThreshold()
        threshold.SetInput(self.surface)
        threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_POINTS", "InteriorLabel")
        threshold.ThresholdBetween(0, 0)
        threshold.Update()
            
        polysurface = vtk.vtkGeometryFilter()
        polysurface.SetInputConnection(threshold.GetOutputPort())
        polysurface.Update()
        self.open_surface = polysurface.GetOutput()

    def get_output(self):
        return self.surface
    
class BoundaryExtractor2d():
    
    
    def __init__(self):
        
        self.surface = None
        self.labels = None
        self.boundary_edges = None
        
    def update(self):
        
        """
        For each marked location on the boudaries extract the start and end point of the 
        boundary edge.
        """
        
        self.boundary_edges = []
        
        for eachLabel in self.labels:
            
            threshold = vtk.vtkThreshold()
            threshold.SetInput(self.surface)
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
                 
            self.boundary_edges.append(region_points)
    
    def get_output(self):
        
        return self.boundary_edges
    
"""
Convert meshes between Tri/Tetgen and Dolfin formats
"""

import numpy as np
import dolfin as df
import vtk

try:
   import cPickle as pickle
except:
   import pickle
    
class DolfinConverter3d():
    
    def __init__(self):
        pass
        
    def generate(self, mesh, regions_set = False):
        editor = df.MeshEditor()
        dolfin_mesh = df.Mesh()
        
        node_locations = mesh.GetNodeLocations()
        connectivity = mesh.GetConnectivity()
        
        editor.open(dolfin_mesh, 3, 3)
        editor.init_vertices(len(node_locations))
        editor.init_cells(len(connectivity))
        
        for i, p in enumerate(node_locations):
            editor.add_vertex(i, np.array(p))
        
        for i, t in enumerate(connectivity):
            editor.add_cell(i, np.array(t, dtype=np.uintp))
        editor.close()
        dolfin_mesh.init()
        return dolfin_mesh
    
    def generate_from_file(self, path, element_label_path, regions_set = False):
        mesh = df.Mesh(path)
        mesh.init() # establish connectivities
        
        with open(element_label_path, 'rb') as handle:
            data  = pickle.load(handle)
        element_labels = data
        
        # set up the vessel and tissue domains
        numcells = len(mesh.cells())
        mesh_func = df.MeshFunction("size_t", mesh, 3)
        mesh_func.set_all(0)
        for idx in range(numcells):
            mesh_func[idx] = int(element_labels[idx])
            
        return mesh, mesh_func
    
class Marker3d():
    
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
    
class FaceBoundary(df.SubDomain):
    
    def set_label(self, label, boundaries, tol = 1.e-3):
        self.label = label
        self.tol = tol
        self.boundary = boundaries[label]
        self.locator = vtk.vtkKdTreePointLocator()
        self.locator.SetDataSet(self.boundary)
        self.points = self.boundary.GetPoints()    
        
    def inside(self, x, on_boundary):
        
        inside = False
        
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
            
        closest_id = self.locator.FindClosestPoint(position)
        dist = np.linalg.norm(np.array(self.points.GetPoint(closest_id)) - position)
        
        return dist <= self.tol
        
def write_meshes(mesh, domain_markers, boundary_markers, output_directory, prefix):
     
    # Save meshes and domains to VTK and XML files
    file = df.File(output_directory + "/" + prefix + "_mesh.pvd")
    file << mesh
      
    file = df.File(output_directory + "/" + prefix + "_domains.pvd")
    file << domain_markers
      
    file = df.File(output_directory + "/" + prefix + "_boundaries.pvd")
    file << boundary_markers
      
    file = df.File(output_directory + "/" + prefix + "_mesh.xml")
    file << mesh
      
    file = df.File(output_directory + "/" + prefix + "_domains.xml")
    file << domain_markers
      
    file = df.File(output_directory + "/" + prefix + "_boundaries.xml")
    file << boundary_markers 