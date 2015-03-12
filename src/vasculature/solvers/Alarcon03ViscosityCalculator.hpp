//
//  Alarcon03ViscosityCalculation.hpp
//  
//
//  Created by Anthony Connor on 01/08/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef _Alarcon03ViscosityCalculation_hpp
#define _Alarcon03ViscosityCalculation_hpp

#include "AbstractVesselPropertyLocalCalculation.hpp"
#include "matrix.hpp"
#include "vectordouble.hpp"
#include <vector>
#include <iostream>
#include "boost/shared_ptr.hpp"

using std::vector;
using std::cout;

template<int DIM, class T>
class Alarcon03ViscosityCalculation : public AbstractVesselPropertyLocalCalculation<DIM,T>
{
    
protected:

    
public:
    
    // constructor
    Alarcon03ViscosityCalculation();

    /**
     *  Virtual destructor.
     */
    virtual ~Alarcon03ViscosityCalculation();
    
    // method for performing the calculation
    virtual void Calculate(boost::shared_ptr<Vessel<DIM,T> > V);
    
};

#endif
