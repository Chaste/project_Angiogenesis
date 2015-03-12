//
//  Alarcon03WallShearStressCalculator.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 07/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include <stdio.h>
#include "Alarcon03WallShearStressCalculator.hpp"

template<int DIM, class T>
Alarcon03WallShearStressCalculator<DIM,T>::Alarcon03WallShearStressCalculator() : AbstractVesselPropertyLocalCalculator<DIM,T>()
{
    
}

template<int DIM, class T>
Alarcon03WallShearStressCalculator<DIM,T>::~Alarcon03WallShearStressCalculator()
{
    
}

template<int DIM, class T>
void Alarcon03WallShearStressCalculator<DIM,T>::Calculate(boost::shared_ptr<Vessel<DIM,T> > V)
{
    
    double WallShearStress;
    
    WallShearStress = 8 * (V->GetViscosity()) * fabs(V->GetFlowRate()) / (PI * pow((V->GetRadius()),3));
    
    V->SetWallShearStress(WallShearStress);
    
}

// Explicit instantiation
template class Alarcon03WallShearStressCalculator<2,int>;
template class Alarcon03WallShearStressCalculator<3,int>;
