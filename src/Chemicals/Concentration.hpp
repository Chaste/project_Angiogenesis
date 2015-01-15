//
//  Concentration.hpp
//  
//
//  Created by Anthony Connor on 28/11/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef ConcentrationDef
#define ConcentrationDef

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using std::string;

class Concentration
{
    
protected:
    
    double Conc;
    string Units;
    
public:
    
    Concentration();
    Concentration(double concentration, string units);
    
    double GetConcentration();
    string GetUnits();
    
    void SetConcentration(double conc,string units);
    void SetConcentration(double conc);
    void SetUnits(string units);
    
};

#endif