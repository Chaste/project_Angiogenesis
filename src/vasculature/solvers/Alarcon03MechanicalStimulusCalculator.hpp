//
//  Alarcon03MechanicalStimulusCalculator.hpp
//  
//
//  Created by Anthony Connor on 01/08/2011.
//  Copyright 2011 Oxford University. All rights reserved.
//

#ifndef _Alarcon03MechanicalStimulusCalculator_hpp
#define _Alarcon03MechanicalStimulusCalculator_hpp

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "CaVascularNetwork.hpp"

template<unsigned DIM>
class Alarcon03MechanicalStimulusCalculator
{
    
protected:
    

	/*
	 * A small constant included to avoid singular behavior at low wall shear stress.
	 */
    double Tau_ref;

    /*
     * The level of wall shear stress expected from the actual intravascular pressure, according
     * to a parametric description of experimental data obtained in the rat mesentry (exhibiting a
     * sigmoidal increase of wall shear stress with increasing pressure.
     */
    double Tau_P;
    
public:
    
    // constructor
    Alarcon03MechanicalStimulusCalculator();
    
    /**
     *  Virtual destructor.
     */
    virtual ~Alarcon03MechanicalStimulusCalculator();
    
    double GetTauP();
    
    double GetTauRef();

    void SetTauRef(double Tau_ref);
    
    // method for performing the Calculation
    /**
        This Calculator has been changed from the original found in Pries1998 in order to better fit experimental data.
        See original paper and relevant test for comparison.
     */
    virtual void Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);
    
    
    
};

#endif
