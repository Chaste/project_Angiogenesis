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
	 * Container for non-spatial node data.
	 */
	boost::shared_ptr<VasculatureData> mpDataContainer;

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
	void AddVessel(boost::shared_ptr<CaVessel<DIM> > vessel);

	/**
       Adds a collection of vessels to the VesselNetwork
	 */
	void AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels);

	/**
       Return true if the input nodes are connected.
	 */
	bool Connected(boost::shared_ptr<VascularNode<DIM> > node1, boost::shared_ptr<VascularNode<DIM> > node2);

	/**
        Return the extents of the vessel network in the form ((xmin, xmax), (ymin, ymax), (zmin, zmax))
	 */
	std::vector<std::pair<double, double> > GetExtents();

	/**
        Return the nodes in the network
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > GetNodes();

	/**
	    Return the only the nodes at the ends of vessels in the network
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > GetVesselEndNodes();

	/**
        Return the vessels in the network
	 */
	std::vector<boost::shared_ptr<CaVessel<DIM> > > GetVessels();

	/**
            Return the vessel at index i in the network
	 */
	boost::shared_ptr<CaVessel<DIM> > GetVessel(unsigned i);

	/**
       Returns the number of vessels contained within the network.
	 */
	unsigned GetNumberOfVessels();

	/**
        Apply the input data to all nodes in the network
	 */
	void SetNodeData(VasculatureData data);

	/**
        Apply the input data to all vessels in the network
	 */
	void SetVesselData(VasculatureData data);

	/*
	 * Translates the network along the provided vector, if a copy is requested the original vessels are copied
	 * (without original non-spatial data) and the new vessels are translated.
	 */
	void Translate(std::vector<double> translation_vector, bool copy = false);

	/**
       Merge nodes with the same spatial location. Useful for
       tidying up networks read from file.
	 */
	void MergeCoincidentNodes();

	/**
     	 Write the VesselNetwork data to a file.
	 */
	void WriteToFile(std::string filename, bool geometry_only = false);

	/**
	 * Returns index of vessel in mVessel.
	 */
	unsigned GetVesselIndex(boost::shared_ptr<CaVessel<DIM> > vessel);

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



	//    /**
	//       Return the Index of the specified vessel
	//     */
	//    unsigned GetVesselIndex(boost::shared_ptr<CaVessel<DIM> > vessel);
	//

	//
	//    /**
	//       Returns the number of nodes contained within the network.
	//     */
	//    unsigned GetNumberOfNodes();
	//
	//    /**
	//       Returns a boost::shared_ptr to this VesselNetwork.
	//
	//       @return boost::shared_ptr<CaVascularNetwork<DIM> >
	//     */
	//    boost::shared_ptr<CaVascularNetwork<DIM> > Shared();
	//
	/**
		Returns a boost::shared_ptr to the node object with the prescribed location in the spatial
		mesh.
	 */
	//boost::shared_ptr<VascularNode<DIM> > GetNode(ChastePoint<DIM> location);

	/**
		Return the number of vessels which occupy a particular location in the spatial mesh.

		@param location ChastePoint<DIM>.
	 */
	//unsigned GetNumberOfVesselsAtLocation(ChastePoint<DIM> location);

	/**
            Returns whether there is a node present at the prescribed location.
	 */
	// bool NodePresentAtLocation(ChastePoint<DIM> location);

	/**
            Returns how many nodes are present at a location.
	 */
	// unsigned NumberOfNodesPresentAtLocation(ChastePoint<DIM> location);

	/**
            Checks whether the prescribed vessel is contained within the vessel network.
	 */
	//bool IsInNetwork(boost::shared_ptr<CaVessel<DIM> > vessel);

	/**
            Checks whether the prescribed node is contained within the vessel network.
	 */
	//bool IsInNetwork(boost::shared_ptr<VascularNode<DIM> > node);

	/**
            Return whether a node is connected to a source node.
	 */
	//bool QueryNodeConnectivity(std::vector<boost::shared_ptr<VascularNode<DIM> > > source_nodes, boost::shared_ptr<VascularNode<DIM> > query_node);

	/**
            Return whether a vector of nodes is connected to a source node.
	 */
	//bool QueryNodeConnectivity(std::vector<boost::shared_ptr<VascularNode<DIM> > > source_nodes, std::vector<boost::shared_ptr<VascularNode<DIM> > > query_nodes);
};

#endif /* CAVASCULARNETWORK_HPP_ */
