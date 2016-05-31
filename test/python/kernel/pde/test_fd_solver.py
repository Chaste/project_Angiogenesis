''' Tests for the hybrid solver module.
'''

import math
import os
import numpy as np
from nose.tools import assert_equals, assert_almost_equals
from casie.utility import assert_almost_equals_lists
import petsc4py, sys
import casie.geometry
import casie.mesh
import casie.pde

        
class TestFiniteDifferenceSolver():
    
    @classmethod
    def setup_class(cls):
        petsc4py.init(sys.argv)
        
    @classmethod
    def teardown_class(self):
        pass
          
    def test_fixed_outer_boundary(self):
        domain = casie.geometry.Part()
        domain.AddCuboid(100.0, 100.0, 100.0)
        
        grid = casie.mesh.RegularGrid()
        grid.GenerateFromPart(domain, 10.0)
        
        pde = casie.pde.HybridLinearEllipticPde()
        pde.SetIsotropicDiffusionConstant(0.003)
        pde.SetContinuumLinearInUTerm(-1.e-5)
        
        bc = casie.pde.HybridBoundaryCondition()
        bc.SetValue(30.0)
        
#         file_handler = casie.core.OutputFileHandler("Casie/TestFiniteDifferenceSolver")
        solver = casie.pde.FiniteDifferenceSolver()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.AddBoundaryCondition(bc)
#         solver.SetFileHandler(file_handler)
#         solver.SetFileName("OuterBoundary")
#         solver.SetWriteSolution(True)
        solver.Solve()
#         print solver.GetVtkSolution()
        
    def test_fixed_left_boundary(self):
        domain = casie.geometry.Part()
        domain.AddCuboid(100.0, 100.0, 10.0)
        
        grid = casie.mesh.RegularGrid()
        grid.GenerateFromPart(domain, 10.0)
        
        pde = casie.pde.HybridLinearEllipticPde()
        pde.SetIsotropicDiffusionConstant(0.003)
        pde.SetContinuumLinearInUTerm(-1.e-5)
        
        for eachFacet in domain.GetFacets():
            if eachFacet.GetCentroid()[0] == 0.0:
                eachFacet.SetData("Boundary", 30.0)
        
        bc = casie.pde.HybridBoundaryCondition()
        bc.SetType(casie.pde.BoundaryConditionType.FACET)
        bc.SetSource(casie.pde.BoundaryConditionSource.LABEL_BASED)
        bc.SetLabelName("Boundary")
        bc.SetDomain(domain)
        
        solver = casie.pde.FiniteDifferenceSolver()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.AddBoundaryCondition(bc)
        solver.Solve()
        
    def test_fixed_left_right_boundary(self):
        domain = casie.geometry.Part()
        domain.AddCuboid(100.0, 100.0, 10.0)
        
        grid = casie.mesh.RegularGrid()
        grid.GenerateFromPart(domain, 10.0)
        
        pde = casie.pde.HybridLinearEllipticPde()
        pde.SetIsotropicDiffusionConstant(0.003)
        pde.SetContinuumLinearInUTerm(-1.e-5)
        
        for eachFacet in domain.GetFacets():
            if eachFacet.GetCentroid()[0] == 0.0:
                eachFacet.SetData("LeftBoundary", 30.0)
        
        for eachFacet in domain.GetFacets():
            if eachFacet.GetCentroid()[0] == 100.0:
                eachFacet.SetData("RightBoundary", 15.0)
        
        left_bc = casie.pde.HybridBoundaryCondition()
        left_bc.SetType(casie.pde.BoundaryConditionType.FACET)
        left_bc.SetSource(casie.pde.BoundaryConditionSource.LABEL_BASED)
        left_bc.SetLabelName("LeftBoundary")
        left_bc.SetDomain(domain)
        
        right_bc = casie.pde.HybridBoundaryCondition()
        right_bc.SetType(casie.pde.BoundaryConditionType.FACET)
        right_bc.SetSource(casie.pde.BoundaryConditionSource.LABEL_BASED)
        right_bc.SetLabelName("RightBoundary")
        right_bc.SetDomain(domain)
        
        solver = casie.pde.FiniteDifferenceSolver()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.AddBoundaryCondition(left_bc)
        solver.AddBoundaryCondition(right_bc)
        solver.Solve()
        