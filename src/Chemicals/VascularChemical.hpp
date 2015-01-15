//
//  VascularChemical.hpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef VascularChemicalDef
#define VascularChemicalDef

#include <iostream>
#include "Concentration.hpp"
#include <cstring>

using std::string;


class VascularChemical
{
    
protected:
    
    string pChemicalName;
    Concentration Conc;
    double PermeabilityOfVesselWallToChemical;
    
public:
    
    VascularChemical(string chemicalname, Concentration conc, double permeabilityofvesselwalltochemical);
    
    string GetChemicalName() const;
    double GetConcentration();
    string GetUnits();
    double GetPermeabilityOfVesselWallToChemical() const;
    
    void SetVascularChemicalConcentration(Concentration conc);
    void SetPermeabilityOfVesselWallToChemical(double permeabilityofvesselwalltochemical);
    
};

#endif