from dolfin import *
import code.settings
import numpy as np
import chaste.geometry.converters
import vtk
    
class OnVesselCentre(SubDomain):
     
    def set_parameters(self, network, dimension, tol = 1.e-3):
        self.network = network
        self.tol = tol
        vtk_to_tri = chaste.geometry.converters.VtkToTri()
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
    
class SolverCg():
    
    def __init__(self, parameters, mesh, domains, boundarys, domain_labels):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.domains = domains
        self.boundaries = boundarys
        self.domain_labels = domain_labels
        
    def run(self):
        
        wall = 3
        
        M = Constant(0.0)
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        
        V = FunctionSpace(self.mesh, 'Lagrange', 1)
        c = TrialFunction(V)
        v = TestFunction(V)
        
        a = inner(nabla_grad(c), nabla_grad(v))*dx - M*c*v*dx
        L = f*v*dx
        
        bc = DirichletBC(V, c_in, self.boundaries, wall)
        
        c = Function(V)
        solve(a == L, c, bc)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(c*dx)
        
        sink_rate = self.parameters["damkohler"] * volume
        
        return c, volume, sink_rate, total_concentration/volume
    
class SolverCgRobin():
    
    def __init__(self, parameters, mesh, domains, boundarys, domain_labels):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.domains = domains
        self.boundaries = boundarys
        self.domain_labels = domain_labels
        
    def run(self):

        wall = 3
        
        M = Constant(0.0)
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        perm = Constant(self.parameters["permeability"])
        
        V = FunctionSpace(self.mesh, 'Lagrange', 1)
        c = TrialFunction(V)
        v = TestFunction(V)
    
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundarys)
    
        a = inner(nabla_grad(c), nabla_grad(v))*dx- M*c*v*dx + perm*c*v*ds(wall)
        L = f*v*dx + perm*c_in*v*ds(wall)
        
        c = Function(V)
        solve(a == L, c)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(c*dx)
        
        sink_rate = self.parameters["damkohler"] * volume
        
        return c, volume, sink_rate, total_concentration/volume

class SolverCgLine():
    
    def __init__(self, parameters, mesh, boundarys, network, dimension):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.network = network
        self.dimension = dimension
        
    def run(self):
        M = Constant(0.0)
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        
        V = FunctionSpace(self.mesh, 'Lagrange', 1)
        c = TrialFunction(V)
        v = TestFunction(V)
    
        a = inner(nabla_grad(c), nabla_grad(v))*dx - M*c*v*dx
        L = f*v*dx
        
        centre = OnVesselCentre()
        centre.set_parameters(self.network, self.dimension)
        bc = DirichletBC(V, c_in, centre, method='pointwise')
        
        c = Function(V)
        solve(a == L, c, bc)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(c*dx)
        
        sink_rate = self.parameters["damkohler"] * volume
        
        return c, volume, sink_rate, total_concentration/volume
    
class SolverCgLineRobin():
    
    def __init__(self, parameters, mesh, boundarys, network, dimension = 3):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.network = network
        self.dimension = dimension
        
    def run(self):
        M = Constant(0.0)
        
        if self.dimension == 3:
            radius = self.network.nodes[0].radius
            S = 2.0 * np.pi * radius
        else:
            S = 2.0
        
        length = 2.0
            
        perm = Constant(self.parameters["permeability"]*S * length)# multiply by vessel surface area
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        
        V = FunctionSpace(self.mesh, 'Lagrange', 1)
        c = TrialFunction(V)
        v = TestFunction(V)
    
        center_domain = VertexFunction("size_t", self.mesh, 0)
        centre = OnVesselCentre()
        centre.set_parameters(self.network, self.dimension)
        centre.mark(center_domain, 1)
        
        dPP = dP(subdomain_data=center_domain)

        a = inner(nabla_grad(c), nabla_grad(v))*dx - M*c*v*dx + perm*c*v*dPP(1)
        L = f*v*dx + perm*c_in*v*dPP(1)
        
        c = Function(V)
        solve(a == L, c)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(c*dx)
        
        sink_rate = self.parameters["damkohler"] * volume
        
        return c, volume, sink_rate, total_concentration/volume