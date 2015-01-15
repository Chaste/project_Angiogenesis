//
//  IntraVascularChemicalCollection.hpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef IntraVascularChemicalCollectionDef
#define IntraVascularChemicalCollectionDef

#include <iostream>
#include <vector>
#include "VascularChemical.hpp"
#include "Concentration.hpp"
#include <cstring>

using std::string;
//using std::vector;

class IntraVascularChemicalCollection
{
    
protected:
    
    std::vector<VascularChemical> ChemCollection;
    
public:
    
    IntraVascularChemicalCollection();
    
    std::vector<VascularChemical>& GetIntraVascularChemicalCollection();
    
    double GetIntraVascularChemicalConcentration(string chemicalname);
    double GetPermeabilityOfVesselWallToChemical(string chemicalname);
    void SetIntraVascularChemicalConcentration(string chemicalname, Concentration conc);
    void SetPermeabilityOfVesselWallToChemical(string chemicalname, double permeabilityofvesselwalltochemical);
    void AddIntraVascularChemical(string chemicalname, Concentration conc, double permeabilityofvesselwalltochemical);

    
};

#endif