import vtk
import dolfin as df
import numpy as np
import chaste.utility.bases as bases

class TriMeshToVtk(bases.SimpleIOBase):
    
    """
    Conversion from a mesh description in the Tri format from Meshpy to
    VTK polydata
    @param self.input a mesh in Tri format from MeshPy
    @return self.output VtkPolydata with element edges as lines
    """
    
    def __init__(self):
        super(TriMeshToVtk, self).__init__()
        
    def update(self):
        points = self.input[0]
        edges = self.input[1]
        self.output = vtk.vtkPolyData() 
        new_vtk_points = vtk.vtkPoints()
        new_vtk_lines = vtk.vtkCellArray()
        
        for eachEdge in edges:
            line = vtk.vtkLine()
            loc1 = (points[eachEdge[0]][0], points[eachEdge[0]][1], 0.0)
            loc2 = (points[eachEdge[1]][0], points[eachEdge[1]][1], 0.0)
            pointId1 = new_vtk_points.InsertNextPoint(loc1)
            pointId2 = new_vtk_points.InsertNextPoint(loc2)
            line.GetPointIds().InsertId(0, pointId1)
            line.GetPointIds().InsertId(1, pointId2)
            new_vtk_lines.InsertNextCell(line)
                
        self.output.SetPoints(new_vtk_points)
        self.output.SetLines(new_vtk_lines)
        return self.output
    
class VtkToTriMesh(bases.SimpleIOBase):
    
    """
    Conversion from vtk polydata lines to MeshPy Tri mesh description
    @param self.input vtkpolydata with lines
    @return self.output a mesh in MeshPy Tri format
    """
    
    def __init__(self):
        super(VtkToTriMesh, self).__init__()
        
    def update(self):
        
        vtk_numPoints = self.input.GetNumberOfPoints()    
        vtk_points = self.input.GetPoints()  
    
        points= []
        edges = []
        for i in range(vtk_numPoints):
            points.append([vtk_points.GetPoint(i)[0], vtk_points.GetPoint(i)[1]])
               
        numCells = self.input.GetNumberOfLines()  
        cellArray = self.input.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(int(seg_id))
            edges.append(point_indices)
        self.output = [points, edges]
        return self.output
    
class ChasteMeshToVtkUnstructured(bases.SimpleIOBase):
    
    """
    Conversion from a mesh in Chaste Tri or Tet format to VTK unstructed grid
    @param self.input a mesh in Chaste Tri or Tet format
    @return self.output VTK unstructed grid
    """
    
    def __init__(self):
        super(ChasteMeshToVtkUnstructured, self).__init__()
        self.dimension = 2
        
    def set_dimension(self, dimension):
        self.dimension = dimension
    
    def update(self):
        self.output = vtk.vtkUnstructuredGrid()
        
        # Add VTK points corresponding to mesh nodes
        points = vtk.vtkPoints()
        locations = self.input.GetNodeLocations()
        points.SetNumberOfPoints(len(locations))
        for idx, eachLocation in enumerate(locations):
            args = [idx] + list(eachLocation)
            points.InsertPoint(*args)
        self.output.SetPoints(points)  

        # Add VTK Tets or Triangles corresponding to mesh elements
        connectivity = self.input.GetConnectivity()
        num_elements = len(connectivity)
        self.output.Allocate(num_elements, num_elements)
        for idx in range(num_elements): 
            if self.dimension == 3:  
                vtkElement = vtk.vtkTetra()
            else:
                vtkElement = vtk.vtkTriangle()

            num_nodes = len(connectivity[idx])
            for jdx in range(num_nodes):  
                vtkElement.GetPointIds().SetId(jdx, connectivity[idx][jdx])                       
            self.output.InsertNextCell(vtkElement.GetCellType(), vtkElement.GetPointIds())
        return self.output
    
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
    
class TriMeshToDolfin(bases.SimpleIOBase):
    
    def __init__(self):
        super(TriMeshToDolfin, self).__init__()
        
        self.boundary_edges = None
        self.domain_edges = None
        self.network = None
        self.surface = None
        self.new_network = None
        self.regions_set = False
        
    def update_centrelines(self, mesh, boundaries):
        
        centres = OnVesselCentre()
        centres.set_network(self.network)
        centres.mark(boundaries, 3)
        
    def update(self):
        editor = df.MeshEditor()
        dolfin_mesh = df.Mesh()
        
        editor.open(dolfin_mesh, 2, 2)
        editor.init_vertices(len(self.input.points))
        editor.init_cells(len(self.input.elements))
        
        element_labels = []
        for i, p in enumerate(self.input.points):
            editor.add_vertex(i, np.array(p))
        
        for i, t in enumerate(self.input.elements):
            editor.add_cell(i, np.array(t, dtype=np.uintp))
            if self.regions_set:
                element_labels.append(self.input.element_attributes[i])
        editor.close()
        dolfin_mesh.init()
        
                # Mark the domains
        mesh_func = df.MeshFunction("size_t", dolfin_mesh, 2)
        mesh_func.set_all(0)
        
        if len(element_labels)>0:
            num_cells = len(dolfin_mesh.cells())
            for idx in range(num_cells):
                mesh_func[idx] = int(element_labels[idx])
    
        # Set up facet markers
        boundaries = df.FacetFunction("size_t", dolfin_mesh)
        boundaries.set_all(0)
        
        if self.network is not None:
            self.update_centrelines(dolfin_mesh, boundaries)
            main_mesh_data = [dolfin_mesh, mesh_func, boundaries]
            
            start_points = []
            for eachRegion in self.boundary_edges:
                for eachEdge in eachRegion:
                    start_points.append((np.array(eachEdge[0]) + np.array(eachEdge[1]))/2.0)
            
            centre_tool = chaste.mesh.centrelines.Extract2d(dolfin_mesh, boundaries, 3)
            centre_tool.seed_points = start_points
            centre_mesh, new_network = centre_tool.generate()
            self.new_network = new_network
            
            r1mesh_func = df.MeshFunction("size_t", centre_mesh, 2)
            r1mesh_func.set_all(0)
            
            r1boundaries = df.FacetFunction("size_t", centre_mesh)
            r1boundaries.set_all(0)
            
            r1_mesh_data = [centre_mesh, r1mesh_func, r1boundaries]
            
            return main_mesh_data, r1_mesh_data, None
        
        facet_count = 0
        two_regions_found = False
        for eachFacet in df.facets(dolfin_mesh):
            if eachFacet.num_entities(2) == 2:
                if mesh_func[int(eachFacet.entities(2)[0])] != mesh_func[int(eachFacet.entities(2)[1])]: 
                    boundaries[facet_count] = 3
                    two_regions_found = True
                elif mesh_func[int(eachFacet.entities(2)[0])] == 1:
                    boundaries[facet_count] = 4
                elif mesh_func[int(eachFacet.entities(2)[0])] == 2: 
                    boundaries[facet_count] = 5
            facet_count += 1
         
        if two_regions_found: 
            r1mesh = df.SubMesh(dolfin_mesh, mesh_func, 1)
            r1mesh_func = df.MeshFunction("size_t", r1mesh, 2)
            r1mesh_func.set_all(0)
            
            r1boundaries = df.FacetFunction("size_t", r1mesh)
            r1boundaries.set_all(0)
            
            r2mesh = df.SubMesh(dolfin_mesh, mesh_func, 2)
            r2mesh_func = df.MeshFunction("size_t", r2mesh, 2)
            r2mesh_func.set_all(0) 
            
            r2boundaries = df.FacetFunction("size_t", r2mesh)
            r2boundaries.set_all(0)
        
        if self.boundary_edges is not None:
            line_marker = LineBoundary()
            line_marker.set_edges(self.boundary_edges[0])
            line_marker.mark(boundaries, 1)
            
            line_marker.set_edges(self.boundary_edges[1])    
            line_marker.mark(boundaries, 2)
            
            line_marker.set_edges(self.boundary_edges[0])
        
            if two_regions_found:
                default_boundary = DefaultBoundary()
                default_boundary.mark(r1boundaries, 3)   
                line_marker.mark(r1boundaries, 1)
                line_marker.set_edges(self.boundary_edges[1])    
                line_marker.mark(boundaries, 2)
                line_marker.mark(r1boundaries, 2)
    
                default_boundary = DefaultBoundary()
                default_boundary.mark(r2boundaries, 3)    
                
                if self.domain_edges is not None:
                    line_marker.set_edges(self.domain_edges)
                    line_marker.mark(r2boundaries, 0)
        
        main_mesh_data = [dolfin_mesh, mesh_func, boundaries]
        if two_regions_found:
            r1_mesh_data = [r1mesh, r1mesh_func, r1boundaries]
            r2_mesh_data = [r2mesh, r2mesh_func, r2boundaries]
        else:
            r1_mesh_data = None
            r2_mesh_data = None   
        self.output = [main_mesh_data, r1_mesh_data, r2_mesh_data]     
        
        return self.output
    
class LineBoundary(df.SubDomain):
    
    def set_edges(self, edges, tol = 1.e-3):
        self.edges = edges
        self.tol = tol
    
    def inside(self, x, on_boundary):
        
        inside = False
        
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
            
        if on_boundary:
            for eachEdge in self.edges:
                if len(eachEdge[0]) == 2:
                    e1loc = np.array((eachEdge[0][0], eachEdge[0][1], 0.0))
                    e2loc = np.array((eachEdge[1][0], eachEdge[1][1], 0.0))
                else:
                    e1loc = np.array(eachEdge[0])
                    e2loc = np.array(eachEdge[1])
                if vtk.vtkLine.DistanceToLine(position, e1loc, e2loc) <= self.tol:
                    dp1 = np.linalg.norm(e1loc - position)
                    dp2 = np.linalg.norm(e2loc - position)
                    dpLine = np.linalg.norm(e1loc - e2loc)
                    
                    if dp1 + dp2 <= dpLine + self.tol:
                        inside = True
                        break
        return inside
    
class DefaultBoundary(df.SubDomain):
    
    def inside(self, x, on_boundary):
        return on_boundary
    
class OnVesselCentre(df.SubDomain):
     
    def set_network(self, network, tol = 1.e-3):
        self.network = network
        self.tol = tol
        vtk_to_tri = chaste.geometry.other.VtkToTri()
        points, edges = vtk_to_tri.generate(self.network)
        self.points = points
        self.edges = edges
        
    def inside(self, x, on_boundary):
         
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
        
        am_inside = False
        for eachEdge in self.edges:
            if len(self.points[eachEdge[0]]) == 2:
                e1loc = np.array((self.points[eachEdge[0]][0], self.points[eachEdge[0]][1], 0.0))
                e2loc = np.array((self.points[eachEdge[1]][0], self.points[eachEdge[1]][1], 0.0))
            else:
                e1loc = np.array(self.points[eachEdge[0]])
                e2loc = np.array(self.points[eachEdge[1]])
            if vtk.vtkLine.DistanceToLine(position, e1loc, e2loc) <= self.tol:
                dp1 = np.linalg.norm(e1loc - position)
                dp2 = np.linalg.norm(e2loc - position)
                dpLine = np.linalg.norm(e1loc - e2loc)
                if dp1 + dp2 <= dpLine + self.tol:
                    am_inside = True
                    break
        return am_inside