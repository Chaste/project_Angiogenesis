#include "RegularGrid.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "SharedPottsMeshGenerator.hpp"
#include "PottsMesh.hpp"

inline int Instantiation()
{
    return
            sizeof(RegularGrid<3, 3>) +
            sizeof(DiscreteContinuumMesh<3, 3>) +
            sizeof(SharedPottsMeshGenerator<3>) +
            sizeof(PottsMesh<3>);
}
