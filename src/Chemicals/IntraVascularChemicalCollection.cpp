//
//  IntraVascularChemicalCollection.cpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#include <iostream>
#include <assert.h>
#include <cstring>
#include "IntraVascularChemicalCollection.hpp"

IntraVascularChemicalCollection::IntraVascularChemicalCollection() : ChemCollection()
{
    
}

std::vector<VascularChemical>& IntraVascularChemicalCollection::GetIntraVascularChemicalCollection()
{
    return ChemCollection;
}

void IntraVascularChemicalCollection::AddIntraVascularChemical(string chemicalname, Concentration conc, double permeabilityofvesselwalltochemical)
{
    
    int addchemical = 1;
    
    for (unsigned i = 0; i < ChemCollection.size(); i++)
    {
        if (ChemCollection[i].GetChemicalName() == chemicalname) 
        {
            SetIntraVascularChemicalConcentration(chemicalname, conc);
            addchemical = 0;
            break;
        }
    }
    
    
    if (addchemical == 1) 
    {
        ChemCollection.push_back(VascularChemical(chemicalname, conc, permeabilityofvesselwalltochemical));
    }
    
    
}



double IntraVascularChemicalCollection::GetIntraVascularChemicalConcentration(string chemicalname)
{
    
    int recognizechemical = 0;
    
    for (unsigned i = 0; i < ChemCollection.size(); i++)
    {
        if (ChemCollection[i].GetChemicalName() == chemicalname)
        {
            recognizechemical = 1;
            return ChemCollection[i].GetConcentration();
        }
    }
    
    if (recognizechemical == 0) 
    {
        std::cout << "Cell does not recognize chemical, " <<  chemicalname  << ".\n";
        assert(recognizechemical == 1);
    }
    
    return 0.0;
    
    
}

double IntraVascularChemicalCollection::GetPermeabilityOfVesselWallToChemical(string chemicalname)
{
    
    int recognizechemical = 0;
    
    for (unsigned i = 0; i < ChemCollection.size(); i++)
    {
        if (ChemCollection[i].GetChemicalName() == chemicalname)
        {
            recognizechemical = 1;
            return ChemCollection[i].GetPermeabilityOfVesselWallToChemical();
        }
    }
    
    if (recognizechemical == 0) 
    {
        std::cout << "Cell does not recognize chemical, " <<  chemicalname  << ".\n";
        assert(recognizechemical == 1);
    }

    return 0.0;
    
}

void IntraVascularChemicalCollection::SetIntraVascularChemicalConcentration(string chemicalname, Concentration conc)
{
    int recognizechemical = 0;
    
    for (unsigned i = 0; i < ChemCollection.size(); i++)
    {
        if (ChemCollection[i].GetChemicalName() == chemicalname)
        {
            
            recognizechemical = 1;
            ChemCollection[i].SetVascularChemicalConcentration(conc);
            
        }
    }
    
    if (recognizechemical == 0) 
    {

        //std::cout << "Vessel does not recognize chemical, " <<  chemicalname   << ". Adding species...\n";
        AddIntraVascularChemical(chemicalname, conc, 0.0);
        //assert(recognizechemical == 1);
    }
    
    
}


void IntraVascularChemicalCollection::SetPermeabilityOfVesselWallToChemical(string chemicalname, double permeabilityofvesselwalltochemical)
{
    int recognizechemical = 0;
    
    for (unsigned i = 0; i < ChemCollection.size(); i++)
    {
        if (ChemCollection[i].GetChemicalName() == chemicalname)
        {
            
            recognizechemical = 1;
            ChemCollection[i].SetPermeabilityOfVesselWallToChemical(permeabilityofvesselwalltochemical);
            
        }
    }
    
    if (recognizechemical == 0) 
    {
        std::cout << "Cell does not recognize chemical, " <<  chemicalname   << "\n";
        assert(recognizechemical == 1);
    }
    
    
}
