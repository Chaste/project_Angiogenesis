//
//  VascularChemical.cpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#include <iostream>
#include <assert.h>
#include <cstring>
#include "VascularChemical.hpp"
#include "Exception.hpp"
#include <sstream>

VascularChemical::VascularChemical(string chemicalname, Concentration conc, double permeabilityofvesselwalltochemical) : pChemicalName(chemicalname),Conc(conc),PermeabilityOfVesselWallToChemical(permeabilityofvesselwalltochemical)
{

}

double VascularChemical::GetPermeabilityOfVesselWallToChemical() const
{
    return PermeabilityOfVesselWallToChemical;    
}


string VascularChemical::GetChemicalName() const
{
    return pChemicalName;
}


double VascularChemical::GetConcentration()
{
    return Conc.GetConcentration();
}


string VascularChemical::GetUnits()
{
    return Conc.GetUnits();
}


void VascularChemical::SetVascularChemicalConcentration(Concentration conc)
{
    try
    {
        if (conc.GetUnits() != Conc.GetUnits())
        {
            std::stringstream sstm;

            sstm.str("");

            sstm << "Vascular chemical concentration units (" << Conc.GetUnits() << ") are not the same as the units of the concentration being prescribed (" << conc.GetUnits() << ").";

            throw Exception(sstm.str(),"VascularChemical.cpp", 57);
        }
        else
        {
            Conc.SetConcentration(conc.GetConcentration(),conc.GetUnits());
        }

    }
    catch (Exception & e)
    {
        std::cout << e.GetMessage() << std::endl;
        e.Terminate("Program exiting","VascularChemical.cpp", 68);
    }

}


void VascularChemical::SetPermeabilityOfVesselWallToChemical(double permeabilityofvesselwalltochemical)
{
    PermeabilityOfVesselWallToChemical = permeabilityofvesselwalltochemical;
}
