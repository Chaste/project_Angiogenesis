//
//  Alarcon03WallShearStressCalculator.hpp
//  
//
//  Created by Anthony Connor on 01/08/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef _Alarcon03WallShearStressCalculator_hpp
#define _Alarcon03WallShearStressCalculator_hpp

#include <vector>
#include <iostream>
#include "boost/shared_ptr.hpp"

#ifndef PI
#define PI 3.14159
#endif

template<int DIM>
class Alarcon03WallShearStressCalculator
{
    
protected:

    
public:
    
    // constructor
    Alarcon03WallShearStressCalculator();
    
    /**
     *  Virtual destructor.
     */
    virtual ~Alarcon03WallShearStressCalculator();
    
    // method for performing the Calculator
    virtual void Calculate(boost::shared_ptr<Vessel<DIM,T> > V);
    


    
};

#endif
