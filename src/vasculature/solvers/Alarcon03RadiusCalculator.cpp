//
//  Alarcon03RadiusCalculation.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 07/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include <stdio.h>
#include "Alarcon03RadiusCalculation.hpp"

template<int DIM,class T>
Alarcon03RadiusCalculation<DIM,T>::Alarcon03RadiusCalculation() : AbstractVesselPropertyLocalCalculation<DIM,T>(),MinRadius(1*pow(10.0,(-6))),MaxRadius(50*pow(10.0,(-6))),TimeStep(0.0001)
{
    
}


template<int DIM,class T>
Alarcon03RadiusCalculation<DIM,T>::~Alarcon03RadiusCalculation()
{
    
}

template<int DIM,class T>
void Alarcon03RadiusCalculation<DIM,T>::SetMinRadius(double min)
{
    MinRadius = min;
    assert(MinRadius >= 1*pow(10.0,(-6)));
}

template<int DIM,class T>
void Alarcon03RadiusCalculation<DIM,T>::SetMaxRadius(double max)
{
    MaxRadius = max;
    assert(MaxRadius >= MinRadius);
    assert(MaxRadius < 0.1);
}

template<int DIM,class T>
void Alarcon03RadiusCalculation<DIM,T>::SetTimestep(double dt)
{
    TimeStep = dt;
}

template<int DIM,class T>
void Alarcon03RadiusCalculation<DIM,T>::Calculate(boost::shared_ptr<Vessel<DIM,T> > V)

{
    
    double TotalStimulus;
    
    TotalStimulus = V->GetMechanicalStimulus() + V->GetMetabolicStimulus() + V->GetUpstreamConductedStimulus() + V->GetDownstreamConductedStimulus() - V->GetShrinkingStimulus();
    
    V->SetRadius(V->GetRadius()*(1.0 + TimeStep*TotalStimulus));
    
    // confine radius to within maximum and minimum vessel radii
    
    if (V->GetRadius() > MaxRadius)
    {
        V->SetRadius(MaxRadius);
    }
    if (V->GetRadius() < MinRadius)
    {
        V->SetRadius(MinRadius);
    }
    
    
}

// Explicit instantiation
template class Alarcon03RadiusCalculation<2,int>;
template class Alarcon03RadiusCalculation<3,int>;
