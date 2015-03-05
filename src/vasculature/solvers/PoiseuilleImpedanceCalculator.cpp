/*
 * PoiseuilleImpedanceCalculator.cpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#include "PoiseuilleImpedanceCalculator.hpp"
#include "CaVesselSegment.hpp"
#include "VasculatureData.hpp"
#include "LinearSystem.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
PoiseuilleImpedanceCalculator<DIM>::PoiseuilleImpedanceCalculator()
{


}

template<unsigned DIM>
PoiseuilleImpedanceCalculator<DIM>::~PoiseuilleImpedanceCalculator()
{


}

template<unsigned DIM>
void PoiseuilleImpedanceCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{


	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();
	std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = vascularNetwork->GetVessels();


	for (unsigned segIndex = 0; segIndex < segments.size(); segIndex++)
	{
		double length = segments[segIndex]->GetLength();
		if (segments[segIndex]->template GetData<double>("Radius") <= 0)
		{
			EXCEPTION("Radius should be a positive number.");
		}
		if (segments[segIndex]->template GetData<double>("Viscosity") <= 0)
		{
			EXCEPTION("Viscosity should be a positive number.");
		}
		double impedance = 8*segments[segIndex]->template GetData<double>("Viscosity")*length/(M_PI*pow(segments[segIndex]->template GetData<double>("Radius"),4.0));

		segments[segIndex]->SetData("Impedance",impedance);
	}

	for (unsigned vesselIndex = 0; vesselIndex < vessels.size(); vesselIndex++)
	{
		double impedance = 0;

		for (unsigned segIndex = 0; segIndex < vessels[vesselIndex]->GetNumberOfSegments(); segIndex++)
		{
			impedance += vessels[vesselIndex]->GetSegment(segIndex)-> template GetData<double>("Impedance");
		}

		vessels[vesselIndex]->SetData("Impedance",impedance);
	}



}

// Explicit instantiation
template class PoiseuilleImpedanceCalculator<2>;
template class PoiseuilleImpedanceCalculator<3>;

