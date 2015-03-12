//
//  Alarcon03ViscosityCalculation.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 07/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include <stdio.h>
#include "Alarcon03ViscosityCalculation.hpp"

template<int DIM,class T>
Alarcon03ViscosityCalculation<DIM,T>::Alarcon03ViscosityCalculation()
: AbstractVesselPropertyLocalCalculation<DIM,T>()
{
    
}

template<int DIM,class T>
Alarcon03ViscosityCalculation<DIM,T>::~Alarcon03ViscosityCalculation()
{
    
}

template<int DIM,class T>
void Alarcon03ViscosityCalculation<DIM,T>::Calculate(boost::shared_ptr<Vessel<DIM,T> > V)
{
    
    double C;
    double MuRel;
    double MuStar45;
    
    double R = pow(10.0,6)*(V->GetRadius()); // scale radius
    double H = V->GetHaematocritLevel();
    double PlasmaViscosity = 3.5*pow(10.0,(-3));
    double Viscosity;
    
    
    C = (0.8 + exp(-0.15*R))*((1.0/(1.0 + pow(10.0,(-11))*pow(2*R,12))) - 1) + (1.0/(1.0 + pow(10.0,(-11))*pow(2*R,12)));
    
    MuStar45 = 6*exp(-0.17*R) + 3.2 - 2.44*exp(-0.06*pow(2*R,0.645));
    
    MuRel = (1 + (MuStar45 - 1)*(((pow((1 - H), C)) - 1)/((pow((1 - 0.45), C)) - 1))*pow((2*R/(2*R - 1.1)),2))*pow((2*R/(2*R - 1.1)),2);
    
    Viscosity = PlasmaViscosity * MuRel;
    
    V->SetViscosity(Viscosity);
    
}

// Explicit instantiation
template class Alarcon03ViscosityCalculation<2,int>;
template class Alarcon03ViscosityCalculation<3,int>;
