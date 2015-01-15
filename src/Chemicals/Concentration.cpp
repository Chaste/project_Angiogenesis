//
//  Concentration.cpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#include <iostream>
#include <assert.h>
#include <cstring>
#include <sstream>
#include "Concentration.hpp"
#include "Exception.hpp"

using std::string;

Concentration::Concentration() : Conc(0.0),Units("")
{
    
}


Concentration::Concentration(double concentration, string units) : Conc(concentration),Units(units)
{
    
}


double Concentration::GetConcentration()
{
    return Conc;
}


string Concentration::GetUnits()
{
    return Units;
}


void Concentration::SetConcentration(double conc,string units)
{

    try
    {
        if (Units != units)
        {
            std::stringstream sstm;

            sstm.str("");

            sstm << "Chemical concentration units (" << Units << ") are not the same as the units of the concentration being prescribed (" << units << ").";

            throw Exception(sstm.str(),"Concentraion.cpp", 54);
        }
        else
        {
            Conc = conc;
        }

    }
    catch (Exception & e)
    {
        std::cout << e.GetMessage() << std::endl;
        e.Terminate("Program exiting","Concentration.cpp", 65);
    }

}

void Concentration::SetConcentration(double conc)
{
    Conc = conc;
}


void Concentration::SetUnits(string units)
{
    Units = units;
}
