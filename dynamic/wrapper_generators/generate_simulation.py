#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):

    builder.class_("OnLatticeSimulationWrapper").include()
    builder.class_("NodeBasedSimulationWrapper").include()
    builder.class_("VascularTumourSolver< 3 >").include()
    builder.class_("VascularTumourModifier< 3 >").include()  
    builder.class_("SimulationManager< 3 >").include()
    
    builder.class_("VascularTumourSolver< 3 >").rename("VascularTumourSolver3")
    builder.class_("VascularTumourModifier< 3 >").rename("VascularTumourModifier3")
    builder.class_("SimulationManager< 3 >").rename("SimulationManager3")
    builder.class_("BetteridgeHaematocritSolver< 3 >").rename("BetteridgeHaematocritSolver3")
    
    return builder
