#include "OnLatticeSimulationWrapper.hpp"
#include "NodeBasedSimulationWrapper.hpp"
#include "VascularTumourSolver.hpp"
#include "VascularTumourModifier.hpp"
#include "SimulationManager.hpp"

template class VascularTumourSolver<3>;
template class VascularTumourModifier<3>;
template class SimulationManager<3>;
