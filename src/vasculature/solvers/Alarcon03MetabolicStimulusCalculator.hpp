//
//  Alarcon03MetabolicStimulusCalculator.hpp
//  
//
//  Created by Anthony Connor on 01/08/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef _Alarcon03MetabolicStimulusCalculator_hpp
#define _Alarcon03MetabolicStimulusCalculator_hpp

#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include "CaVascularNetwork.hpp"

template<unsigned DIM>
class Alarcon03MetabolicStimulusCalculator
{
    
protected:

    double Q_ref;
    double k_m;
    double MaxStimulus;
    
public:
    
    // constructor
    Alarcon03MetabolicStimulusCalculator();
    
    /**
     *  Virtual destructor.
     */
    virtual ~Alarcon03MetabolicStimulusCalculator();

    double GetQRef();

    double GetKm();

    double GetMaxStimulus();

    void SetQRef(double qref);
    
    void SetKm(double km);
    
    void SetMaxStimulus(double maxstimulus);
    
    // method for performing the calculation
    virtual void Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);

};

#endif
