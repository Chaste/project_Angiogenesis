import numpy as np
from dolfin import *
import chaste.population.vessel
import chaste.simulation
import vtk
try:
   import cPickle as pickle
except:
   import pickle

class vel_field_3d(Expression):
    def eval(self, value, x):
        value[0] = 0.0
        value[1] = 0.0
        value[2] = 200.0 * 2.0 * (1.0 - (x[0]*x[0] + x[1]*x[1])/(10.0*10.0))
    def value_shape(self):
        return (3,)
    
class vel_field_2d(Expression):
    
    def __init__(self, degree, edges, mesh, tol):
        
        self.edges = edges
        self.mesh = mesh
        self.tol = tol
    
    def eval_cell(self, value, x, ufc_cell):
        
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
        
        vel = 0.0
        for eachEdge in self.edges:
            
            pd_edge1 = np.array((eachEdge[0][0], eachEdge[0][1], 0.0))
            pd_edge2 = np.array((eachEdge[1][0], eachEdge[1][1], 0.0))
            if vtk.vtkLine.DistanceToLine(position, pd_edge1, pd_edge2) <= self.tol:
                dp1 = np.linalg.norm(pd_edge1 - position)
                dp2 = np.linalg.norm(pd_edge2 - position)
                length = np.linalg.norm(pd_edge1 - pd_edge2)
                
                if dp1 + dp2 <= length + self.tol:
                    mid = (pd_edge1 + pd_edge2) /2.0
                    dist = np.linalg.norm(mid - position)
                    vel = 50.0*2.0*(1.0 - dist*dist/((length/2.0)*(length/2.0)))
                
        value[0] = -vel*n[0]
        value[1] = -vel*n[1]
        
    def value_shape(self):
        return (2,)
    
class VesselInlet(SubDomain):
     
    def set_parameters(self, network, dimension, tol = 1.e-3):
        self.network = network
        self.dimension = dimension#
        self.tol = tol
 
    def inside(self, x, on_boundary):
        im_inside = x[1] < 1.e-6
        return im_inside

class Solver():
    
    def __init__(self, parameters, mesh, boundarys, dimension = 3, edges=None):
        self.mesh = mesh
        self.parameters = parameters
        self.boundaries = boundarys
        self.dimension = dimension
        self.edges = edges
        
    def run(self):
    
        wall = 3
        inlet = 1
        
        v_in = self.parameters["average_velocity"]
        mu = Constant(self.parameters["effective_viscosity"])
    
            # Define function spaces
        P2 = VectorElement("Lagrange", self.mesh.ufl_cell(), 2)
        P1 = FiniteElement("Lagrange", self.mesh.ufl_cell(), 1)
        TH = P2 * P1
        W = FunctionSpace(self.mesh, TH)
        
        # No Slip Boundary Condition on vessel wall
        if self.dimension == 2:
            no_slip = Constant((0,0))
        else:
            no_slip = Constant((0,0,0))
            
        bc0 = DirichletBC(W.sub(0), no_slip, self.boundaries, wall)
            
        # Positive flow in y direction at the inlet
        if self.dimension == 2:
            inflow = vel_field_2d(degree=2, mesh = self.mesh, edges = self.edges, tol=1.e-3)
        else:
            inflow = vel_field_3d(degree=2, mesh = self.mesh, edges = self.edges, tol=1.e-3)
            
        bc1 = DirichletBC(W.sub(0), inflow, self.boundaries, inlet)
    
        bcs = [bc0, bc1]
        
        # Specify the problem
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        
        if self.dimension == 2:        
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))
        a = (mu * inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
        L = inner(f, v)*dx
        
        # Compute solution
        w = Function(W)
        solve(a == L, w, bcs)
        
        # Split the mixed solution using deepcopy
        u, p = w.split(True)
        
        return [u, p]
    
class Solver1d():
    
    def __init__(self, parameters, mesh, dimension = 3):
        self.mesh = mesh
        self.parameters = parameters
        self.dimension = dimension
        
    def run(self):

        v_in = self.parameters["average_velocity"]
        mu = Constant(self.parameters["effective_viscosity"])
    
            # Define function spaces
        P2 = VectorFunctionSpace(self.mesh, "Lagrange", 2)
        P1 = FunctionSpace(self.mesh, "Lagrange", 1)
        W = P2 * P1
        
        # No Slip Boundary Condition on vessel wall
        #if self.dimension == 2:
        #    no_slip = Constant((0,0))
        #else:
        #    no_slip = Constant((0,0,0))
            
        #bc0 = DirichletBC(W.sub(0), no_slip, self.boundaries, wall)
        
        # Positive flow in y direction at the inlet
        if self.dimension == 2:
            inflow = Constant((0.0, v_in))
        else:
            inflow = Constant((0.0, 0.0, v_in))
            
        #inflow = vel_field(degree=2)
        
        inlet = VesselInlet()
        #inlet.set_parameters(self.network, self.dimension)
        bc1 = DirichletBC(W.sub(0), inflow, inlet)
        
        #bc1 = DirichletBC(W.sub(0), inflow, self.boundaries, inlet)
    
        bcs = [bc1]
        
        # Specify the problem
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        
        if self.dimension == 2:        
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))
        a = (mu * inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
        L = inner(f, v)*dx
        
        # Compute solution
        w = Function(W)
        #solve(a == L, w, bcs, solver_parameters={"linear_solver": "mumps"})
        solve(a == L, w, bcs)#, solver_parameters={"linear_solver": "umfpack"})
        
        # Split the mixed solution using deepcopy
        u, p = w.split(True)
        
        return [u, p]