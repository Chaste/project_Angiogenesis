/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "CaVessel.hpp"

template<unsigned DIM>
CaVesselSegment<DIM>::CaVesselSegment(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2)
: mNodes(std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > (pNode1, pNode2)),
  mDataContainer(),
  mId(0),
  mLabel(""),
  mVessel(boost::weak_ptr<CaVessel<DIM> >()),
  mRadius(1.0),
  mHaematocrit(0.45),
  mFlowRate(0.0),
  mImpedance(0.0)
  {

  }

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> >CaVesselSegment<DIM>::Create(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2)
{
	boost::shared_ptr<CaVesselSegment<DIM> > pSelf(new CaVesselSegment<DIM>(pNode1, pNode2));

	// Make sure the specified nodes are not the same
	if (pNode1 == pNode2)
	{
		EXCEPTION("Attempted to assign the same node to both ends of a vessel segment.");
	}

	// Add the segment to the nodes
	pNode1->AddSegment(pSelf->shared_from_this());
	pNode2->AddSegment(pSelf->shared_from_this());
	return pSelf;
}

template<unsigned DIM>
CaVesselSegment<DIM>::~CaVesselSegment()
{
}

template<unsigned DIM>
void CaVesselSegment<DIM>::AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel)
{
	mVessel = pVessel;
}

template<unsigned DIM>
const VasculatureData& CaVesselSegment<DIM>::rGetDataContainer() const
{
	return mDataContainer;
}

template<unsigned DIM>
template<typename T> T CaVesselSegment<DIM>::GetData(const std::string& variableName)
{
	return mDataContainer.GetData<T>(variableName);
}

template<unsigned DIM>
std::vector<std::string> CaVesselSegment<DIM>::GetDataKeys(bool castable_to_double) const
{
	return mDataContainer.GetKeys(castable_to_double);
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetDistance(const ChastePoint<DIM>& rPoint)
{
	ChastePoint<DIM> location1 = GetNode(0)->GetLocation();
	ChastePoint<DIM> location2 = GetNode(1)->GetLocation();

	std::vector<double> segment_vector(DIM);
	std::vector<double> point_vector(DIM);

	// Get a vector along the segment and one from the point to one end of the segment
	for(unsigned i=0; i<DIM; i++)
	{
		segment_vector[i] = location2[i] - location1[i];
		point_vector[i] = rPoint[i] - location1[i];
	}

	// Get segment_vector.point_vector
	double dot_product_segment_point = 0.0;
	for(unsigned i=0; i<DIM; i++)
	{
		dot_product_segment_point += segment_vector[i] * point_vector[i];
	}

	if(dot_product_segment_point <= 0.0)
	{
		return GetNode(0)->GetDistance(rPoint);
	}

	// Get point_vector.point_vector
	double dot_product_segment_segment = 0.0;
	for(unsigned i=0; i<DIM; i++)
	{
		dot_product_segment_segment += segment_vector[i] * segment_vector[i];
	}

	if(dot_product_segment_segment <= dot_product_segment_point)
	{
		return GetNode(1)->GetDistance(rPoint);
	}

	double projection_ratio = dot_product_segment_point / dot_product_segment_segment;
	std::vector<double> projected_point(DIM);

	for(unsigned i=0; i<DIM; i++)
	{
		projected_point[i] += location1[i] + projection_ratio * segment_vector[i];
	}

	double distance_squared;

	distance_squared = pow(projected_point[0]-rPoint[0], 2) + pow(projected_point[1]-rPoint[1], 2);
	if(DIM==3)
	{
		distance_squared += pow(projected_point[2]-rPoint[2],2);
	}

	return std::sqrt(distance_squared);
}

template<unsigned DIM>
unsigned CaVesselSegment<DIM>::GetId() const
{
	return mId;
}

template<unsigned DIM>
const std::string& CaVesselSegment<DIM>::rGetLabel() const
{
	return mLabel;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetLength() const
{
	ChastePoint<DIM> point1 = mNodes.first->GetLocation();
	ChastePoint<DIM> point2 = mNodes.second->GetLocation();

	double length_squared = pow(point2[0] - point1[0], 2) + pow(point2[1] - point1[1], 2);

	if (DIM == 3u)
	{
		length_squared += pow(point2[2] - point1[2], 2);
	}
	return std::sqrt(length_squared);
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetRadius() const
{
	return mRadius;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetHaematocrit() const
{
	return mHaematocrit;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetFlowRate() const
{
	return mFlowRate;
}

template<unsigned DIM>
double CaVesselSegment<DIM>::GetImpedance() const
{
	return mImpedance;
}

template<unsigned DIM>
ChastePoint<DIM> CaVesselSegment<DIM>::GetMidPoint()
{
	ChastePoint<DIM> point1 = mNodes.first->GetLocation();
	ChastePoint<DIM> point2 = mNodes.second->GetLocation();

	if (DIM ==2)
	{
		ChastePoint<DIM> mid_point((point2[0] + point1[0])/2.0, (point2[1] + point1[1])/2.0);
		return mid_point;
	}
	else
	{
		ChastePoint<DIM> mid_point((point2[0] + point1[0])/2.0, (point2[1] + point1[1])/2.0, (point2[2] + point1[2])/2.0);
		return mid_point;
	}
}

template<unsigned DIM>
boost::shared_ptr<VascularNode<DIM> > CaVesselSegment<DIM>::GetNode(unsigned index)
{
	if (index == 0u)
	{
		return mNodes.first;
	}
	else if (index == 1u)
	{
		return mNodes.second;
	}
	else
	{
		EXCEPTION("A node index other than 0 or 1 has been requested for a Vessel Segment.");
	}
}

template<unsigned DIM>
std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > CaVesselSegment<DIM>::GetNodes()
{
	return mNodes;
}

template<unsigned DIM>
ChastePoint<DIM> CaVesselSegment<DIM>::GetPointProjection(const ChastePoint<DIM>& rPoint)
{
	ChastePoint<DIM> location1 = GetNode(0)->GetLocation();
	ChastePoint<DIM> location2 = GetNode(1)->GetLocation();

	std::vector<double> segment_vector(DIM);
	std::vector<double> point_vector(DIM);

	// Get a vector along the segment and one from the point to one end of the segment
	for(unsigned i=0; i<DIM; i++)
	{
		segment_vector[i] = location2[i] - location1[i];
		point_vector[i] = rPoint[i] - location1[i];
	}

	// Get segment_vector.point_vector
	double dot_product_segment_point = 0.0;
	for(unsigned i=0; i<DIM; i++)
	{
		dot_product_segment_point += segment_vector[i] * point_vector[i];
	}

	// Get point_vector.point_vector
	double dot_product_segment_segment = 0.0;
	for(unsigned i=0; i<DIM; i++)
	{
		dot_product_segment_segment += segment_vector[i] * segment_vector[i];
	}

	if(dot_product_segment_point <= 0.0 || dot_product_segment_segment <= dot_product_segment_point)
	{
		EXCEPTION("Projection of point is outside segment.");
	}

	double projection_ratio = dot_product_segment_point / dot_product_segment_segment;
	std::vector<double> projected_point(DIM);

	for(unsigned i=0; i<DIM; i++)
	{
		projected_point[i] += location1[i] + projection_ratio * segment_vector[i];
	}

	return ChastePoint<DIM>(projected_point);
}

template<unsigned DIM>
ChastePoint<DIM> CaVesselSegment<DIM>::GetUnitTangent()
{

	ChastePoint<DIM> location1 = GetNode(0)->GetLocation();
	ChastePoint<DIM> location2 = GetNode(1)->GetLocation();

	std::vector<double> tangent;
	double length_squared = 0.0;
	for(unsigned i=0; i<DIM; i++)
	{
		double distance = location2[i] - location1[i];
		length_squared += distance * distance;
		tangent.push_back(distance);
	}

	double length = std::sqrt(length_squared);
	if (length <= 1.e-12)
	{
		EXCEPTION("Segment is too short to reliably calculate a unit tangent.");
	}
	for(unsigned i=0; i<DIM; i++)
	{
		tangent[i] /= length;
	}

	return ChastePoint<DIM>(tangent);
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVesselSegment<DIM>::GetVessel()
{
	if(mVessel.lock())
	{
		return mVessel.lock();
	}
	else
	{
		EXCEPTION("A vessel has been requested but this segment doesn't have one.");
	}
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::HasDataKey(const std::string& rKey) const
{
	return mDataContainer.HasKey(rKey);
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::HasNode(boost::shared_ptr<VascularNode<DIM> > pNode)
{
	return(pNode == GetNode(0) || pNode == GetNode(1));
}

template<unsigned DIM>
bool CaVesselSegment<DIM>::IsConnectedTo(boost::shared_ptr<CaVesselSegment<DIM> > otherSegment)
{

	bool isConnectedToSegment = false;

	if (otherSegment == Shared())
	{
		EXCEPTION("Vessel segment cannot be connected to itself.");
	}

	if (this->GetNode(0) == otherSegment->GetNode(0) && this->GetNode(1) == otherSegment->GetNode(1))
	{
		EXCEPTION("Vessel segments should not have identical nodes at both ends (and thus be completely overlapping).");
	}
	if (this->GetNode(1) == otherSegment->GetNode(0) && this->GetNode(0) == otherSegment->GetNode(1))
	{
		EXCEPTION("Vessel segments should not have identical nodes at both ends (and thus be completely overlapping).");
	}

	if (this->GetNode(0) == otherSegment->GetNode(0)
			|| this->GetNode(0) == otherSegment->GetNode(1)
			|| this->GetNode(1) == otherSegment->GetNode(0)
			|| this->GetNode(1) == otherSegment->GetNode(1))
	{
		isConnectedToSegment = true;
	}

	return isConnectedToSegment;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::RemoveVessel()
{
	mVessel = boost::weak_ptr<CaVessel<DIM> >();
}

template<unsigned DIM>
void CaVesselSegment<DIM>::ReplaceNode(unsigned oldNodeIndex, boost::shared_ptr<VascularNode<DIM> >  pNewNode)
{
	if (oldNodeIndex == 0u)
	{
		mNodes.first->RemoveSegment(Shared());
		mNodes.first = pNewNode;
		mNodes.first->AddSegment(Shared());
	}
	else if (oldNodeIndex == 1u)
	{
		mNodes.second->RemoveSegment(Shared());
		mNodes.second = pNewNode;
		mNodes.second->AddSegment(Shared());
	}
	else
	{
		EXCEPTION("A node index other than 0 or 1 has been requested for a Vessel Segment.");
	}
}

template<unsigned DIM>
template<typename T> void CaVesselSegment<DIM>::SetData(const std::string& variableName, T value)
{
	mDataContainer.SetData(variableName, value);
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetDataContainer(const VasculatureData& rDataContainer)
{
	mDataContainer = rDataContainer;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetId(unsigned id)
{
	mId = id;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetLabel(const std::string& rLabel)
{
	mLabel = rLabel;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetRadius(double radius)
{
	mRadius = radius;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetHaematocrit(double haematocrit)
{
	mHaematocrit = haematocrit;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetFlowRate(double flowRate)
{
	mFlowRate = flowRate;
}

template<unsigned DIM>
void CaVesselSegment<DIM>::SetImpedance(double impedance)
{
	mImpedance = impedance;
}

template<unsigned DIM>
boost::shared_ptr<CaVesselSegment<DIM> > CaVesselSegment<DIM>::Shared()
{
	boost::shared_ptr<CaVesselSegment<DIM> > pSegment = this->shared_from_this();
	return pSegment;
}

// Explicit instantiation
template class CaVesselSegment<2>;
template class CaVesselSegment<3>;

template bool CaVesselSegment<2>::GetData<bool>(const std::string& variableName);
template double CaVesselSegment<2>::GetData<double>(const std::string& variableName);
template unsigned CaVesselSegment<2>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> CaVesselSegment<2>::GetData<std::vector<double> >(const std::string& variableName);
template void CaVesselSegment<2>::SetData(const std::string& variableName, bool value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, double value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, unsigned value);
template void CaVesselSegment<2>::SetData(const std::string& variableName, std::vector<double> value);

template bool CaVesselSegment<3>::GetData<bool>(const std::string& variableName);
template double CaVesselSegment<3>::GetData<double>(const std::string& variableName);
template unsigned CaVesselSegment<3>::GetData<unsigned>(const std::string& variableName);
template std::vector<double> CaVesselSegment<3>::GetData<std::vector<double> >(const std::string& variableName);
template void CaVesselSegment<3>::SetData(const std::string& variableName, bool value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, double value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, unsigned value);
template void CaVesselSegment<3>::SetData(const std::string& variableName, std::vector<double> value);
