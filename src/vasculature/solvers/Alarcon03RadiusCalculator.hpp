//
//  Alarcon03RadiusCalculation.hpp
//  
//
//  Created by Anthony Connor on 01/08/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef _Alarcon03RadiusCalculation_hpp
#define _Alarcon03RadiusCalculation_hpp

#include "AbstractVesselPropertyLocalCalculation.hpp"
#include "matrix.hpp"
#include "vectordouble.hpp"
#include <vector>
#include <iostream>
#include "boost/shared_ptr.hpp"

using std::vector;
using std::cout;

template<int DIM, class T>
class Alarcon03RadiusCalculation : public AbstractVesselPropertyLocalCalculation<DIM,T>
{
    
protected:

    double MinRadius;
    double MaxRadius;
    double TimeStep;
    
public:
    
    // constructor
    Alarcon03RadiusCalculation();
    
    /**
     *  Virtual destructor.
     */
    virtual ~Alarcon03RadiusCalculation();
    
    void SetMinRadius(double min);
    
    void SetMaxRadius(double max);

    void SetTimestep(double dt);
    
    // method for performing the calculation
    virtual void Calculate(boost::shared_ptr<Vessel<DIM,T> > V);
    
};

#endif
