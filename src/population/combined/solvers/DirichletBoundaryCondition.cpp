/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "Facet.hpp"
#include "DirichletBoundaryCondition.hpp"
//#include "SimpleCell.hpp"

template<unsigned DIM>
DirichletBoundaryCondition<DIM>::DirichletBoundaryCondition()
    :   mpNetwork(),
//        mpCellPopulation(),
        mpDomain(),
        mPoints(),
        mType(BoundaryConditionType::POINT),
        mSource(BoundaryConditionSource::PRESCRIBED),
        mLabel("Default"),
        mValue(0.0)
{

}

template<unsigned DIM>
DirichletBoundaryCondition<DIM>::~DirichletBoundaryCondition()
{

}

template<unsigned DIM>
boost::shared_ptr<DirichletBoundaryCondition<DIM> > DirichletBoundaryCondition<DIM>::Create()
{
    MAKE_PTR(DirichletBoundaryCondition<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
BoundaryConditionType::Value DirichletBoundaryCondition<DIM>::GetType()
{
    return mType;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

//template<unsigned DIM>
//void DirichletBoundaryCondition<DIM>::SetCellPopulation(boost::shared_ptr<SimpleCellPopulation<DIM> > pCellPopulation)
//{
//    mpCellPopulation = pCellPopulation;
//}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpDomain = pDomain;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetPoint(c_vector<double, DIM> point)
{
    mPoints.push_back(point);
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetPoints(std::vector<c_vector<double, DIM> > points)
{
    mPoints.insert(mPoints.end(), points.begin(), points.end());
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetType(BoundaryConditionType::Value boundaryType)
{
    mType = boundaryType;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetSource(BoundaryConditionSource::Value boundarySource)
{
    mSource = boundarySource;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetLabelName(const std::string& label)
{
    mLabel = label;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetValue(double value)
{
    mValue = value;
}

template<unsigned DIM>
std::pair<bool, double> DirichletBoundaryCondition<DIM>::GetValue(c_vector<double, DIM> location, double tolerance)
{
    // Check the boundary condition type
    if(mType == BoundaryConditionType::OUTER)
    {
        std::pair<bool, double>(true, mValue);
    }

    else if(mType == BoundaryConditionType::POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of boundary condition");
        }
        else
        {
            if(norm_2(location-mPoints[0]) < tolerance)
            {
                return std::pair<bool, double>(true, mValue);
            }
            else
            {
                return std::pair<bool, double>(false, mValue);
            }
        }
    }
    else if(mType == BoundaryConditionType::MULTI_POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of boundary condition");
        }
        else
        {
            for(unsigned idx=0; idx<mPoints.size(); idx++)
            {
                if(norm_2(location-mPoints[0]) < tolerance)
                {
                    return std::pair<bool, double>(true, mValue);
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == BoundaryConditionType::FACET)
    {
        if(!mpDomain)
        {
            EXCEPTION("A part is required for this type of boundary condition");
        }
        else
        {
            std::vector<boost::shared_ptr<Facet> > facets =  mpDomain->GetFacets();
            for(unsigned idx=0; idx<facets.size();idx++)
            {
                if(facets[idx]->ContainsPoint(location))
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, facets[idx]->GetData(mLabel));
                    }
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == BoundaryConditionType::VESSEL_LINE)
    {
        if(!mpNetwork)
        {
            EXCEPTION("A vessel network is required for this type of boundary condition");
        }
        else
        {
            std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
            for (unsigned idx = 0; idx <  segments.size(); idx++)
            {
                if (segments[idx]->GetDistance(location) <= tolerance)
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, segments[idx]->template GetData<double>(mLabel));
                    }
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == BoundaryConditionType::VESSEL_VOLUME)
    {
        if(!mpNetwork)
        {
            EXCEPTION("A vessel network is required for this type of boundary condition");
        }
        else
        {
            std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
            for (unsigned idx = 0; idx <  segments.size(); idx++)
            {
                if (segments[idx]->GetDistance(location) <= segments[idx]->GetRadius() + tolerance)
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, segments[idx]->template GetData<double>(mLabel));
                    }
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == BoundaryConditionType::CELL_POINT)
    {
//        if(!mpCellPopulation)
//        {
//            EXCEPTION("A simple cell population is required for this type of boundary condition");
//        }
//
//        std::vector<boost::shared_ptr<SimpleCell<DIM> > > cells = mpCellPopulation->GetCells();
//        for(unsigned idx=0; idx<cells.size(); idx++)
//        {
//            if (norm_2(cells[idx]->rGetLocation()-location)<tolerance)
//            {
//                return std::pair<bool, double>(true, mValue);
//            }
//        }
        return std::pair<bool, double>(false, mValue);
    }
    else if(mType == BoundaryConditionType::IN_PART)
    {
        if(!mpDomain)
        {
            EXCEPTION("A part is required for this type of boundary condition");
        }
        else
        {

            if(mpDomain->IsPointInPart(location))
            {
                if(BoundaryConditionSource::PRESCRIBED)
                {
                    return std::pair<bool, double>(true, mValue);
                }
                else
                {
                    return std::pair<bool, double>(true, mValue);
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    return std::pair<bool, double>(false, mValue);
}

template<unsigned DIM>
std::vector<std::pair<bool, double> > DirichletBoundaryCondition<DIM>::GetValues(std::vector<c_vector<double, DIM> > locations, double tolerance)
{
    std::vector<std::pair<bool, double> > result;
    for(unsigned idx=0; idx<locations.size();idx++)
    {
        result.push_back(GetValue(locations[idx],tolerance));
    }
    return result;
}

// Explicit instantiation
template class DirichletBoundaryCondition<2>;
template class DirichletBoundaryCondition<3>;
