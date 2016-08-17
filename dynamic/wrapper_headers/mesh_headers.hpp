#include "RegularGrid.hpp"
#include "HybridMesh.hpp"
#include "SharedPottsMeshGenerator.hpp"
#include "PottsMesh.hpp"

inline int Instantiation()
{
    return
            sizeof(RegularGrid<3, 3>) +
            sizeof(HybridMesh<3, 3>) +
            sizeof(SharedPottsMeshGenerator<3>) +
            sizeof(PottsMesh<3>);
}
