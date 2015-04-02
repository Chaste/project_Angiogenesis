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

#ifndef CAVASCULARNETWORK_HPP_
#define CAVASCULARNETWORK_HPP_

#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <map>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the boost deprecated warning for now (gcc4.3)
#include <boost/config.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#ifdef CHASTE_VTK
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK
#include "CaVessel.hpp"
#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "VasculatureData.hpp"
#include "GeometryTransform.hpp"

/*
    A vessel network is a collection of vessels
 */

template<unsigned DIM>
class CaVascularNetwork : public boost::enable_shared_from_this<CaVascularNetwork<DIM> >
{

private:

	/**
       Container for Vessels in the VesselNetwork.
	 */
	std::vector<boost::shared_ptr<CaVessel<DIM> > > mVessels;

	/**
	 * Container for non-spatial vessel network data.
	 */
	VasculatureData mDataContainer;

public:

	/*
	 * Constructor
	 */
	CaVascularNetwork();

	/*
	 * Destructor
	 */
	~CaVascularNetwork();

	/**
        Adds a vessel to the VesselNetwork.
	 */
	void AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel);

	/**
       Adds a collection of vessels to the VesselNetwork
	 */
	void AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels);

	/**
        Get the node nearest to the specified location
	 */
	boost::shared_ptr<VascularNode<DIM> > GetNearestNode(ChastePoint<DIM>& rLocation);

	/**
         Get the segment nearest to the specified location
	 */
	boost::shared_ptr<CaVesselSegment<DIM> > GetNearestSegment(const ChastePoint<DIM>& rLocation);

	/**
         Get the segment nearest to the specified location
	 */
	boost::shared_ptr<CaVessel<DIM> > GetNearestVessel(const ChastePoint<DIM>& rLocation);

	/**
        Return the extents of the vessel network in the form ((xmin, xmax), (ymin, ymax), (zmin, zmax))
	 */
	std::vector<std::pair<double, double> > GetExtents();

	/**
        Return the nodes in the network
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > GetNodes();

	/**
        Return the node with the prescribed index
	 */
	boost::shared_ptr<VascularNode<DIM> > GetNode(unsigned index);

	/**
        Return the node with the prescribed index
	 */
	unsigned GetNodeIndex(boost::shared_ptr<VascularNode<DIM> > node);

	/**
        Return the nodes in the network in the form of a vector
	 */
	std::vector<boost::shared_ptr<VascularNode<DIM> > > GetVectorOfNodes();

	/**
	   Returns the number of nodes in the network.
	 */
	unsigned GetNumberOfNodes();

	/**
	   Returns the number of vessels in the network.
	 */
	unsigned GetNumberOfVessels();

	/**
        Return the vessel with the specified index in the network
	 */
	boost::shared_ptr<CaVessel<DIM> > GetVessel(unsigned index);

	/**
	    Return the only the nodes at the ends of vessels in the network
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > GetVesselEndNodes();

	/**
	   Return the Index of the specified vessel
	 */
	unsigned GetVesselIndex(boost::shared_ptr<CaVessel<DIM> > pVessel);

	/**
       Return the Index of the specified vessel segment
	 */
	unsigned GetVesselSegmentIndex(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

	/**
       Return the vessel segments in the network
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetVesselSegments();

	/**
        Return the vessels in the network
	 */
	std::vector<boost::shared_ptr<CaVessel<DIM> > > GetVessels();

	/**
        Return whether a node is connected to a source node.
	 */
	bool IsConnected(boost::shared_ptr<VascularNode<DIM> > pSourceNode, boost::shared_ptr<VascularNode<DIM> > pQueryNode);

	/**
        Return whether a vector of nodes is connected to a vector of source nodes.
	 */
	std::vector<bool> IsConnected(std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes, std::vector<boost::shared_ptr<VascularNode<DIM> > > queryNodes);

	/**
	 * Return whether node is in network.
	 */
	bool NodeIsInNetwork(boost::shared_ptr<VascularNode<DIM> > pSourceNode);

	/**
       Merge nodes with the same spatial location. Useful for
       tidying up networks read from file.
	 */
	void MergeCoincidentNodes();

	/**
        Apply the input data to all nodes in the network
	 */
	void SetNodeData(VasculatureData data);

	/**
        Apply the input data to all vessel segments in the network
	 */
	void SetVesselData(VasculatureData data);

	/**
	        Apply the input data to all vessels in the network
	 */
	void SetSegmentData(VasculatureData data);

	/*
	 * Translates the network along the provided vector, if a copy is requested the original vessels are copied
	 * (without original non-spatial data) and the new vessels are translated.
	 */
	void Translate(const std::vector<double>& rTranslationVector, bool copy = false);

	/**
     	 Write the VesselNetwork data to a file.
	 */
	void Write(const std::string& rFilename, bool geometryOnly = false);

	/**
	 * Outputs connectivity of vessels to file in graphviz format (.gv).
	 */
	void VisualiseVesselConnectivity(std::string output_filename);

private:

	/**
        Create a deep copy of all existing vessels and return a set of newly added nodes.
        This is useful for performing geometric transformations on the network.
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > DeepCopyVessels();

};

#endif /* CAVASCULARNETWORK_HPP_ */
