#include "DiscreteSource.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "DiscreteContinuumNonLinearEllipticPde.hpp"
#include "AbstractDiscreteContinuumSolver.hpp"
#include "AbstractRegularGridDiscreteContinuumSolver.hpp"
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
            return  sizeof(AbstractDiscreteContinuumSolver<3>) +
                    sizeof(AbstractRegularGridDiscreteContinuumSolver<3>) +
                    sizeof(FiniteDifferenceSolver<3>) +
                    sizeof(FiniteElementSolver<3>) +
                    sizeof(DistanceMap<3>) +
                    sizeof(FunctionMap<3>) +
                    sizeof(GreensFunctionSolver<3>) +
                    sizeof(CellStateDependentDiscreteSource<3>) +
                    sizeof(DiscreteContinuumBoundaryCondition<3>) +
                    sizeof(DiscreteContinuumNonLinearEllipticPde<3, 3>) +
                    sizeof(DiscreteContinuumLinearEllipticPde<3, 3>) +
                    sizeof(AbstractLinearEllipticPde<3, 3>) +
                    sizeof(DiscreteSource<3>);
        }
    }
}
