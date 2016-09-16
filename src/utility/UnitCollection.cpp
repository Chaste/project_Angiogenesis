#include "UnitCollection.hpp"

template<class T>
UnitTester<T>::UnitTester() :
    mMyMass()
{
}

template<class T>
UnitTester<T>::~UnitTester()
{

}

//void UnitTester::SetMass(units::quantity<unit::mass> inputMass)
//{
//    mMyMass =  3.0 * inputMass;
//}

template<class T>
units::quantity<T> UnitTester<T>::GetMass()
{
    return mMyMass;
}

//explicit inst.
template class UnitTester<unit::mass>;

