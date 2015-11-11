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

template<unsigned DIM>
DirichletBoundaryCondition<DIM>::DirichletBoundaryCondition()
    :   mpNetwork(),
        mpCellPopulation(),
        mpDomain(),
        mPoints(),
        mType(BoundaryConditionType::POINT),
        mSource(BoundaryConditionSource::PRESCRIBED),
        mLabel("Default"),
        mValue(0.0),
        mpRegularGrid(),
        mpMesh()
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
double DirichletBoundaryCondition<DIM>::GetValue()
{
    return mValue;
}

template<unsigned DIM>
BoundaryConditionType::Value DirichletBoundaryCondition<DIM>::GetType()
{
    return mType;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::UpdateBoundaryConditionContainer(boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> > pContainer)
{
    double node_distance_tolerance = 1.e-3;

    bool apply_boundary = true;
    bool use_boundry_nodes = false;

    if(mType == BoundaryConditionType::OUTER)
    {
        pContainer->DefineConstantDirichletOnMeshBoundary(mpMesh.get(), mValue);
        apply_boundary = false;
    }
    else if(mType == BoundaryConditionType::FACET || mType == BoundaryConditionType::VESSEL_VOLUME)
    {
        use_boundry_nodes = true;
    }

    if(apply_boundary)
    {
        if(!use_boundry_nodes)
        {
            typename TetrahedralMesh<DIM, DIM>::NodeIterator iter = mpMesh->GetNodeIteratorBegin();
            while (iter != mpMesh->GetNodeIteratorEnd())
            {
                std::pair<bool,double> result = GetValue((*iter).GetPoint().rGetLocation(), node_distance_tolerance);
                if(result.first)
                {
                    ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(result.second);
                    pContainer->AddDirichletBoundaryCondition(&(*iter), p_fixed_boundary_condition, 0, false);
                }
                ++iter;
            }
        }
        else
        {
            typename TetrahedralMesh<DIM, DIM>::BoundaryNodeIterator iter = mpMesh->GetBoundaryNodeIteratorBegin();
            while (iter < mpMesh->GetBoundaryNodeIteratorEnd())
            {
                std::pair<bool,double> result = GetValue((*iter)->GetPoint().rGetLocation(), node_distance_tolerance);
                if(result.first)
                {
                    ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(result.second);
                    pContainer->AddDirichletBoundaryCondition(*iter, p_fixed_boundary_condition);
                }
                ++iter;
            }
        }
    }
}

template<unsigned DIM>
std::pair<bool, double> DirichletBoundaryCondition<DIM>::GetValue(c_vector<double,DIM> location, double tolerance)
{
    std::pair<bool, double> result(false, 0.0);
    if(mType == BoundaryConditionType::POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of boundary condition");
        }
        else
        {
            for(unsigned jdx=0; jdx<mPoints.size(); jdx++)
            {
                if(norm_2(location-mPoints[jdx]) < tolerance)
                {
                    return std::pair<bool, double>(true, mValue);
                }
            }
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
            for(unsigned jdx=0; jdx<facets.size();jdx++)
            {
                if(facets[jdx]->ContainsPoint(location) && (facets[jdx]->GetData("Boundary")>0.0))
                {
                    return std::pair<bool, double>(true, mValue);
                }
            }
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
            for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
            {
                if (segments[jdx]->GetDistance(location) <= tolerance)
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, segments[jdx]->template GetData<double>(mLabel));
                    }
                }
            }
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
            for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
            {
                if (segments[jdx]->GetDistance(location) <= segments[jdx]->GetRadius() + tolerance)
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, segments[jdx]->template GetData<double>(mLabel));
                    }
                }
            }
        }
    }

    else if(mType == BoundaryConditionType::CELL)
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
                    return  std::pair<bool, double>(true, mValue);
                }
                else
                {
                    return  std::pair<bool, double>(true, mValue);
                }
            }
        }
    }
    return result;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::UpdateRegularGridBoundaryConditions(boost::shared_ptr<std::vector<std::pair<bool, double> > >pBoundaryConditions,
                                                                          double tolerance)
{
    if(! mpRegularGrid)
    {
        EXCEPTION("A grid has not been set for the determination of boundary condition values. For FE solvers use GetBoundaryConditionContainer()");
    }

    // Check the boundary condition type
    if(mType == BoundaryConditionType::OUTER)
    {
        for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
        {
            (*pBoundaryConditions)[idx] = std::pair<bool, double> (mpRegularGrid->IsOnBoundary(idx), mValue);
        }
    }
    else if(mType == BoundaryConditionType::POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of boundary condition");
        }
        else
        {
            for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
            {
                for(unsigned jdx=0; jdx<mPoints.size(); jdx++)
                {
                    if(norm_2(mpRegularGrid->GetLocationOf1dIndex(idx)-mPoints[jdx]) < tolerance)
                    {
                        (*pBoundaryConditions)[idx] = std::pair<bool, double>(true, mValue);
                        break;
                    }
                }
            }
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
            for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
            {
                std::vector<boost::shared_ptr<Facet> > facets =  mpDomain->GetFacets();
                for(unsigned jdx=0; jdx<facets.size();jdx++)
                {
                    if(facets[jdx]->ContainsPoint(mpRegularGrid->GetLocationOf1dIndex(idx)))
                    {
                        if(BoundaryConditionSource::PRESCRIBED)
                        {
                            (*pBoundaryConditions)[idx] = std::pair<bool, double>(true, mValue);
                            break;
                        }
                        else
                        {
                            (*pBoundaryConditions)[idx] = std::pair<bool, double>(true, facets[jdx]->GetData(mLabel));
                            break;
                        }
                    }
                }
            }
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
            for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
            {
                std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
                for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
                {
                    if (segments[jdx]->GetDistance(mpRegularGrid->GetLocationOf1dIndex(idx)) <= tolerance)
                    {
                        if(BoundaryConditionSource::PRESCRIBED)
                        {
                            (*pBoundaryConditions)[idx] = std::pair<bool, double>(true, mValue);
                            break;
                        }
                        else
                        {
                            (*pBoundaryConditions)[idx] = std::pair<bool, double>(true, segments[jdx]->template GetData<double>(mLabel));
                            break;
                        }
                    }
                }
            }
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
            for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
            {
                std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
                for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
                {
                    if (segments[jdx]->GetDistance(mpRegularGrid->GetLocationOf1dIndex(idx)) <= segments[jdx]->GetRadius() + tolerance)
                    {
                        if(BoundaryConditionSource::PRESCRIBED)
                        {
                            (*pBoundaryConditions)[idx] =  std::pair<bool, double>(true, mValue);
                            break;
                        }
                        else
                        {
                            (*pBoundaryConditions)[idx] =  std::pair<bool, double>(true, segments[jdx]->template GetData<double>(mLabel));
                            break;
                        }
                    }
                }
            }
        }
    }

    else if(mType == BoundaryConditionType::CELL)
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
    }
    else if(mType == BoundaryConditionType::IN_PART)
    {
        if(!mpDomain)
        {
            EXCEPTION("A part is required for this type of boundary condition");
        }
        else
        {
            for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
            {
                if(mpDomain->IsPointInPart(mpRegularGrid->GetLocationOf1dIndex(idx)))
                {
                    if(BoundaryConditionSource::PRESCRIBED)
                    {
                        (*pBoundaryConditions)[idx] =  std::pair<bool, double>(true, mValue);
                        break;
                    }
                    else
                    {
                        (*pBoundaryConditions)[idx] =  std::pair<bool, double>(true, mValue);
                        break;
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation)
{
    mpCellPopulation = boost::shared_ptr<AbstractCellPopulation<DIM> >(&rCellPopulation);
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpDomain = pDomain;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetPoints(std::vector<c_vector<double, DIM> > points)
{
    mPoints = points;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetSource(BoundaryConditionSource::Value boundarySource)
{
    mSource = boundarySource;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetType(BoundaryConditionType::Value boundaryType)
{
    mType = boundaryType;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pRegularGrid)
{
    mpRegularGrid = pRegularGrid;
}

template<unsigned DIM>
void DirichletBoundaryCondition<DIM>::SetMesh(boost::shared_ptr<TetrahedralMesh<DIM, DIM> > pMesh)
{
    mpMesh = pMesh;
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
void DirichletBoundaryCondition<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

// Explicit instantiation
template class DirichletBoundaryCondition<2>;
template class DirichletBoundaryCondition<3>;
