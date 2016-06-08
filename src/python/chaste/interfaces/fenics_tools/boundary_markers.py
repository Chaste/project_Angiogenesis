import numpy as np
import vtk
import dolfin as df

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