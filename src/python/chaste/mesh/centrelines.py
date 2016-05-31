from dolfin import *
import vtk
import numpy as np
import casie.population.vessel
    
class Extract2d():
    
    def __init__(self, mesh, boundary_markers, boundary_label):
        self.mesh = mesh
        self.boundary_markers = boundary_markers
        self.boundary_label = boundary_label
        self.seed_points = None
        self.new_network = None
        self.num_divisions = 5
        
    def generate(self):
        
        # Get the centreline mesh
        facet_count = 0
        centre_facets = []
        unique_verts = []
        vert_coords = []
        
        for eachFacet in facets(self.mesh):
            if self.boundary_markers[facet_count] == self.boundary_label:
                v0 = Vertex(self.mesh, eachFacet.entities(0)[0])
                v1 = Vertex(self.mesh, eachFacet.entities(0)[1])
                centre_facets.append((v0, v1))
                
                if v0.index() not in unique_verts:
                    unique_verts.append(v0.index())
                    vert_coords.append(v0)
                if v1.index() not in unique_verts:
                    unique_verts.append(v1.index())
                    vert_coords.append(v1)   
            facet_count += 1      
            
        # construct the points 
        editor = MeshEditor()
        centre_mesh = Mesh()
        
        editor.open(centre_mesh, 1, 2)
        editor.init_vertices(len(unique_verts) + self.num_divisions*len(self.seed_points))
        editor.init_cells(len(centre_facets) + self.num_divisions*len(self.seed_points))
        
        polyData = vtk.vtkPolyData() 
        new_vtk_points = vtk.vtkPoints()
        new_vtk_lines = vtk.vtkCellArray()

        for idx, eachPoint in enumerate(vert_coords):
            editor.add_vertex(idx, np.array((eachPoint.x(0), eachPoint.x(1)), dtype=np.float))
            new_vtk_points.InsertNextPoint((eachPoint.x(0), eachPoint.x(1), 0.0))
            
        for idx, eachFacet in enumerate(centre_facets):
            id1 = unique_verts.index(eachFacet[0].index())
            id2 = unique_verts.index(eachFacet[1].index())
            editor.add_cell(idx, np.array((id1, id2), dtype=np.uintp))
            line = vtk.vtkLine()
            line.GetPointIds().InsertId(0, id1)
            line.GetPointIds().InsertId(1, id2)
            new_vtk_lines.InsertNextCell(line)

        polyData.SetPoints(new_vtk_points)
        polyData.SetLines(new_vtk_lines)
        
        vtk_numPoints = polyData.GetNumberOfPoints()  
        # Add extensions
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(polyData)
        
        point_counter = len(vert_coords)
        element_counter = len(centre_facets)
        num_points = 5
        for idx, eachLoc in enumerate(self.seed_points):
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            
            end_loc = np.array(eachLoc)
            mesh_loc = np.array((new_vtk_points.GetPoint(closest_id)[0], new_vtk_points.GetPoint(closest_id)[1]))
            vector = mesh_loc - end_loc
            length = np.linalg.norm(vector)
            norm_vector = vector/length
            spacing = length/float(num_points)
            for jdx in range(num_points):
                loc = end_loc + norm_vector*spacing*float(jdx)
                new_vtk_points.InsertNextPoint(probe_loc)
                editor.add_vertex(point_counter, loc)
                if jdx>0:
                    new_vtk_lines.InsertNextCell(2)
                    new_vtk_lines.InsertCellPoint(point_counter)
                    new_vtk_lines.InsertCellPoint(point_counter-1)
                    editor.add_cell(element_counter, np.array((point_counter, point_counter-1), dtype=np.uintp))
                    element_counter+=1
                point_counter+=1
            new_vtk_lines.InsertNextCell(2)
            new_vtk_lines.InsertCellPoint(point_counter-1)
            new_vtk_lines.InsertCellPoint(closest_id)
            editor.add_cell(element_counter, np.array((point_counter-1, closest_id), dtype=np.uintp))
            element_counter+=1
            
        editor.close()
        centre_mesh.init()
        centre_mesh.order()
        
        new_polyData = vtk.vtkPolyData() 
        new_polyData.SetPoints(new_vtk_points)
        new_polyData.SetLines(new_vtk_lines)
        
        return centre_mesh, new_polyData