#include "UnitCollection.hpp"

UnitTester::UnitTester() :
    mMyMass(10.0*unit::kg)
{
}

UnitTester::~UnitTester()
{

}

void UnitTester::SetMass(units::quantity<unit::mass> inputMass)
{
    mMyMass =  3.0 * inputMass;
}

units::quantity<unit::mass> UnitTester::GetMass()
{
    return mMyMass;
}

