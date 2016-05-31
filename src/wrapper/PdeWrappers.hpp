#ifndef PDEWRAPPERS_HPP_
#define PDEWRAPPERS_HPP_

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include "SmartPointers.hpp"
#include "FunctionMap.hpp"
#include "AbstractRegularGridHybridSolver.hpp"
#include "AbstractHybridSolver.hpp"
#include "RegularGrid.hpp"

using namespace boost::python;

typedef AbstractHybridSolver<3> AbstractHybridSolver3;
struct AbstractHybridSolverWrap : AbstractHybridSolver3, wrapper<AbstractHybridSolver3>
{
    void Solve()
    {
        this->get_override("Solve")();
    }

    void Update()
    {
        this->get_override("Update")();
    }

    void Setup()
    {
        this->get_override("Setup")();
    }

    void Write()
    {
        this->get_override("Write")();
    }

    std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, 3> > samplePoints)
    {
        this->get_override("GetSolutionAtPoints")();
    }

    std::vector<double> GetSolutionAtGridPoints(boost::shared_ptr<RegularGrid<3> > pGrid)
    {
        this->get_override("GetSolutionAtGridPoints")();
    }

    void UpdateCellData()
    {
        this->get_override("UpdateCellData")();
    }
};

typedef AbstractRegularGridHybridSolver<3> AbstractRegularGridHybridSolver3;
struct AbstractRegularGridHybridSolverWrap : AbstractRegularGridHybridSolver3, wrapper<AbstractRegularGridHybridSolver3>
{
    void Solve()
    {
        this->get_override("Solve")();
    }

    void Update()
    {
        this->get_override("Update")();
    }

    vtkSmartPointer<vtkImageData> GetVtkSolution()
    {
        if (override GetVtkSolution = this->get_override("GetVtkSolution"))
        {
            return GetVtkSolution();
        }
        return AbstractRegularGridHybridSolver<3>::GetVtkSolution();
    }

    vtkSmartPointer<vtkImageData> default_GetVtkSolution()
    {
        return this->AbstractRegularGridHybridSolver<3>::GetVtkSolution();
    }
};

#endif /* PDEWRAPPERS_HPP_ */
