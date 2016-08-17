#include "DiscreteSource.hpp"
#include "HybridBoundaryCondition.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridNonLinearEllipticPde.hpp"
#include "AbstractHybridSolver.hpp"
#include "AbstractRegularGridHybridSolver.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "FiniteElementSolver.hpp"
#include "DistanceMap.hpp"
#include "FunctionMap.hpp"
#include "GreensFunctionSolver.hpp"
#include "CellStateDependentDiscreteSource.hpp"

namespace chaste
{
    namespace pde
    {
        inline int Instantiation()
        {
            return  sizeof(AbstractHybridSolver<3>) +
                    sizeof(AbstractRegularGridHybridSolver<3>) +
                    sizeof(FiniteDifferenceSolver<3>) +
                    sizeof(FiniteElementSolver<3>) +
                    sizeof(DistanceMap<3>) +
                    sizeof(FunctionMap<3>) +
                    sizeof(GreensFunctionSolver<3>) +
                    sizeof(CellStateDependentDiscreteSource<3>) +
                    sizeof(HybridBoundaryCondition<3>) +
                    sizeof(HybridNonLinearEllipticPde<3, 3>) +
                    sizeof(HybridLinearEllipticPde<3, 3>) +
                    sizeof(AbstractLinearEllipticPde<3, 3>) +
                    sizeof(DiscreteSource<3>);
        }
    }
}
