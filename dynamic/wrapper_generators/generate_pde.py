#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):

    builder.class_("DiscreteSource< 3 >").include()
    builder.class_("HybridBoundaryCondition< 3 >").include()
    builder.class_("HybridNonLinearEllipticPde< 3, 3 >").include()
    builder.class_("HybridLinearEllipticPde< 3, 3 >").include()
    builder.class_("AbstractLinearEllipticPde< 3, 3 >").include()
    builder.class_("AbstractHybridSolver< 3 >").include()
    builder.class_("AbstractRegularGridHybridSolver< 3 >").include()
    builder.class_("FiniteDifferenceSolver< 3 >").include()
    builder.class_("FiniteElementSolver< 3 >").include()
    builder.class_("FunctionMap< 3 >").include()
    builder.class_("GreensFunctionSolver< 3 >").include()
    builder.class_("CellStateDependentDiscreteSource< 3 >").include()
    
    builder.class_("DiscreteSource< 3 >").rename("DiscreteSource3")
    builder.class_("HybridNonLinearEllipticPde< 3, 3 >").rename("HybridNonLinearEllipticPde3")
    builder.class_("HybridLinearEllipticPde< 3, 3 >").rename("HybridLinearEllipticPde3")
    builder.class_("AbstractLinearEllipticPde< 3, 3 >").rename("AbstractLinearEllipticPde3")
    builder.class_("AbstractHybridSolver< 3 >").rename("AbstractHybridSolver3")
    builder.class_("AbstractRegularGridHybridSolver< 3 >").rename("AbstractRegularGridHybridSolver3")
    builder.class_("FiniteDifferenceSolver< 3 >").rename("FiniteDifferenceSolver3")
    builder.class_("FiniteElementSolver< 3 >").rename("FiniteElementSolver3")
    builder.class_("DistanceMap< 3 >").rename("DistanceMap3")
    builder.class_("FunctionMap< 3 >").rename("FunctionMap3")
    builder.class_("GreensFunctionSolver< 3 >").rename("GreensFunctionSolver3")
    builder.class_("CellStateDependentDiscreteSource< 3 >").rename("CellStateDependentDiscreteSource3")
    
    builder.class_("BoundaryConditionType").include()
    builder.class_("BoundaryConditionSource").include()
    builder.class_("SourceType").include()
    builder.class_("SourceStrength").include()  

    return builder