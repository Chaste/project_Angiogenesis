#include "VesselNetworkReader.hpp"
#include "NodeFlowProperties.hpp"
#include "SegmentFlowProperties.hpp"
#include "VesselFlowProperties.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "AbstractVesselNetworkComponent.hpp"
#include "VesselNetworkWriter.hpp"

template class NodeFlowProperties<3>;
template class SegmentFlowProperties<3>;
template class VesselFlowProperties<3>;
template class VesselNode<3>;
template class VesselSegment<3>;
template class Vessel<3>;
template class VesselNetwork<3>;
template class VesselNetworkWriter<3>;
template class AbstractVesselNetworkComponent<3>;
template class VasculatureGenerator<3>;
template class VesselNetworkReader<3>;
