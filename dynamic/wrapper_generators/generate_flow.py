#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):
    
    builder.class_("FlowSolver< 3 >").include()
    builder.class_("WallShearStressCalculator< 3 >").include()
    builder.class_("VesselImpedanceCalculator< 3 >").include()
    builder.class_("BetteridgeHaematocritSolver< 3 >").include()
    
    builder.class_("FlowSolver< 3 >").rename("FlowSolver3")
    builder.class_("WallShearStressCalculator< 3 >").rename("WallShearStressCalculator3")
    builder.class_("VesselImpedanceCalculator< 3 >").rename("VesselImpedanceCalculator3")
    builder.class_("BetteridgeHaematocritSolver< 3 >").rename("BetteridgeHaematocritSolver3")
    
    return builder
