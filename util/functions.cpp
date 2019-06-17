/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-02-17
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "functions.h"

#include <string>
#include <algorithm>    // std::transform
#include <stack>
#include <utility>
#include <fstream>
#include <limits>

#include "lib_disc/domain.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_grid/grid/grid_base_object_traits.h"
#include "lib_grid/global_attachments.h" // global attachments
#include "lib_grid/tools/surface_view.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/algorithms/element_side_util.h" // GetOpposingSide
#include "lib_disc/domain_util.h" // LoadDomain
#include "pcl/pcl_base.h"


namespace ug {
namespace cable_neuron {



template <typename TDomain>
void scale_domain(SmartPtr<TDomain> dom, number scale)
{
	// get attachment accessors
	typename TDomain::position_accessor_type& aaPos = dom->position_accessor();
	ANumber aDiam(GlobalAttachments::attachment<ANumber>("diameter"));
	UG_COND_THROW(!dom->grid()->has_vertex_attachment(aDiam),
		"Diameter attachment not available.");
	Grid::VertexAttachmentAccessor<ANumber> aaDiam(*dom->grid(), aDiam);

	UG_COND_THROW(dom->grid()->num_levels() > 1,
		"Grid must not be refined for this function to work properly.")

	// iterate vertices
	VertexIterator it = dom->grid()->template begin<Vertex>(0);
	VertexIterator itEnd = dom->grid()->template end<Vertex>(0);
	for (; it != itEnd; ++it)
	{
		Vertex* vrt = *it;

		aaPos[vrt] *= scale;
		aaDiam[vrt] *= scale;
	}
}




enum ErrorConstants
{
	EC_NO_ERROR = 0,
	EC_CYCLIC = 1,
	EC_DOUBLE_INDEX = 1 << 1,
	EC_MISSING_INDEX = 1 << 2
};
typedef Flag<ErrorConstants, unsigned int, EC_NO_ERROR>	ErrorState;




template <typename TDomain>
bool is_acyclic(SmartPtr<TDomain> dom, int verbosity)
{
	if (verbosity > 1)
		UG_LOGN("Checking for cycles in the geometry.")

#ifdef UG_PARALLEL
	// TODO: Enhance! (This should be an error, but only if any neuron truly is cut by partitioning.)
	if (pcl::NumProcs() > 1)
		UG_LOGN("The function checking whether the domain is acyclic is intended to be used in serial mode "
				"and may not give correct results otherwise.\n"
				"Please run it on a single processor or make sure that no neuron is cut by domain partitioning.");
#endif

	typedef typename Grid::traits<Edge>::secure_container edge_list;
	typename TDomain::position_accessor_type& aaPos = dom->position_accessor();

	EdgeIterator iter = dom->grid()->template begin<Edge>(0);
	EdgeIterator iterEnd = dom->grid()->template end<Edge>(0);
	std::stack<std::pair<Edge*, Vertex*> > stack;

	// marker for treated vertices
	dom->grid()->begin_marking();

	// init stack with the two edges connected to the first vertex
	for (; iter != iterEnd; ++iter)
	{
		// look for edges disconnected from all previous edges
		if (dom->grid()->is_marked(*iter)) continue;

		// push both ends of initial edge to stack
		stack.push(std::make_pair(*iter, (*iter)->vertex(0)));
		stack.push(std::make_pair(*iter, (*iter)->vertex(1)));

		while (!stack.empty())
		{
			Edge* e = stack.top().first;
			Vertex* v = stack.top().second;
			stack.pop();

			// mark edge
			dom->grid()->mark(e);

			// add any (other) edge connected to vertex
			edge_list el;
			dom->grid()->associated_elements(el, v);
			for (size_t k = 0; k < el.size(); ++k)
			{
				if (el[k] != e)
				{
					// return false if already marked (as this is proof of a cycle)
					if (dom->grid()->is_marked(el[k]))
					{
						if (verbosity > 1)
							UG_LOGN("Found cycle involving vertex at " << aaPos[v] << ".");

						dom->grid()->end_marking();
						return false;
					}

					Vertex* o;
					el[k]->get_opposing_side(v, &o);

					stack.push(std::make_pair(el[k], o));
				}
			}
		}
	}

	dom->grid()->end_marking();

	return true;
}



template <typename TDomain>
int check_presyn_indices(SmartPtr<TDomain> dom, int verbosity)
{
	if (verbosity > 1)
		UG_LOGN("Checking presynaptic indices.")

	ErrorState err_state;

	UG_COND_THROW(!dom.valid(), "Domain invalid.\n");
	SmartPtr<MultiGrid> grid = dom->grid();
	UG_COND_THROW(!grid.valid(), "There is no grid associated to the domain.");

	// handle presynapse index attachment
	if (!GlobalAttachments::is_declared("presyn_index"))
		UG_THROW("GlobalAttachment 'presyn_index' not available.");

	AUInt aPSI = GlobalAttachments::attachment<AUInt>("presyn_index");
	if (!grid->has_vertex_attachment(aPSI)) grid->attach_to_vertices(aPSI);

	Grid::VertexAttachmentAccessor<AUInt> aaPSI = Grid::VertexAttachmentAccessor<AUInt>(*grid, aPSI);

// presynaptic vertex indexing
	typedef typename TDomain::subset_handler_type sh_type;
	SmartPtr<sh_type> sh = dom->subset_handler();
	UG_COND_THROW(!sh.valid(), "The subset handler assigned to the domain taken from the passed\n"
				  "CableEquation object is not valid.")

	// find number of presynaptic indices
	// on each processor, get the maximal index, then communicate and resize.
	SurfaceView sv(sh);
	GridLevel gl;
	typedef SurfaceView::traits<Vertex>::const_iterator it_type;
	it_type it = sv.begin<Vertex>(gl, SurfaceView::MG_ALL);
	it_type it_end = sv.end<Vertex>(gl, SurfaceView::MG_ALL);

	uint max = 0;
	std::list<uint> psi_list;
	for (; it != it_end; ++it)
	{
		Vertex* vrt = *it;
		uint idx = aaPSI[vrt];
		if (idx != (uint) -1)
		{
			psi_list.push_back(idx);
			if (idx > max) max = idx;
		}
	}

	// communicate
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		unsigned long localMax = max;
		com.allreduce(&localMax, &max, 1, PCL_DT_UNSIGNED_LONG, PCL_RO_MAX);
	}
#endif

	std::vector<size_t> psi_present(max+1, 0);

	// now check whether all indices are present
	std::list<uint>::iterator psi_it = psi_list.begin();
	std::list<uint>::iterator psi_it_end = psi_list.end();

	for (; psi_it != psi_it_end; ++psi_it)
		++psi_present[*psi_it];

	// communicate (all-to-rank0)
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1 && psi_present.size())
	{
		pcl::ProcessCommunicator com;

		size_t* localData;
		size_t sz = psi_present.size();
		localData = new size_t[sz];
		memcpy(localData, &psi_present[0], sizeof(size_t)*sz);
		com.reduce(localData, &psi_present[0],
					  sz, pcl::template DataTypeTraits<size_t>::get_data_type(), PCL_RO_SUM, 0);
		delete[] localData;
	}

	if (pcl::ProcRank() == 0)
	{
#endif

	for (size_t i = 0; i < max+1; ++i)
	{
		if (!psi_present[i])
		{
			err_state |= EC_MISSING_INDEX;
			if (verbosity > 1)
				UG_LOGN("Missing presynaptic index " << i << ".")
		}
		else if (psi_present[i] > 1)
		{
			err_state |= EC_DOUBLE_INDEX;
			if (verbosity > 1)
				UG_LOGN("Presynaptic index " << i << "is not unique.")
		}
	}

#ifdef UG_PARALLEL
	}
#endif

	return err_state.get();
}


template <typename TDomain>
int check_domain(SmartPtr<TDomain> dom, int verbosity)
{
	ErrorState err_state;

	try
	{
		// check for cyclicity
		if (!is_acyclic(dom, verbosity))
			err_state |= EC_CYCLIC;

		// check presynapse index attachments
		err_state |= check_presyn_indices(dom, verbosity);
	}
	UG_CATCH_THROW("An error occurred within one of the tests in the domain.");

	// output errors
	if (verbosity > 0)
	{
		UG_LOGN("");
		UG_LOGN("-----------------------------------------------------");
		UG_LOGN("-- Domain check results                            --");
		UG_LOGN("--                                                 --");
		if (err_state.contains(EC_CYCLIC))
			UG_LOGN("-- Domain contains a cycle.                        --");
		if (err_state.contains(EC_DOUBLE_INDEX))
			UG_LOGN("-- At least one presynaptic index is not unique.   --");
		if (err_state.contains(EC_MISSING_INDEX))
			UG_LOGN("-- At least one presynaptic index is missing.      --");
		if (err_state.get() == EC_NO_ERROR)
			UG_LOGN("-- Everything checks out.                          --")

		UG_LOGN("-----------------------------------------------------");
		UG_LOGN("");
	}

	return err_state.get();
}


// This is a debugging tool for a very specific problem.
// Do not compile but for the specific debuging purpose.
#if 0
template <typename TDomain>
void test_vertices(SmartPtr<TDomain> dom)
{
	UG_COND_THROW(!dom.valid(), "Domain invalid.\n");
	SmartPtr<MultiGrid> grid = dom->grid();
	UG_COND_THROW(!grid.valid(), "There is no grid associated to the domain.");

	typename TDomain::position_accessor_type& aaPos = dom->position_accessor();

	VertexIterator it = grid->template begin<Vertex>(0);
	size_t vptr_val = *reinterpret_cast<size_t*>(*it);

	VertexIterator it_end = grid->template end<Vertex>(0);
	size_t cnt = 1;
	for (++it; it != it_end; ++it)
	{
		if (*reinterpret_cast<size_t*>(*it) != vptr_val)
		{
			UG_THROW("Vertex with deviating virtual table pointer:\n"
					"Vertex " << cnt << "/" << grid->num_vertices()
					<< " at " << aaPos[*it] << ".\n"
					"Correct vPtr is " << std::hex << std::setfill('0') << std::setw(16)
					<< vptr_val << ", deviating vPtr is " << std::setfill('0') << std::setw(16) << ".");
		}
		++cnt;
	}

	UG_LOGN("Correct vPtr is " << std::hex << std::setfill('0')
			<< std::setw(16) << vptr_val << std::dec);
}
#endif


static
void depth_first_search
(
	Grid& g,
	Grid::VertexAttachmentAccessor<Attachment<uint> >& aaNID,
	Vertex* v,
	uint id
)
{
	if (aaNID[v] != (uint) -1) return;
    aaNID[v] = id;

    Grid::traits<Edge>::secure_container edges;
    g.associated_elements(edges, v);

    size_t sz = edges.size();
    for (size_t i = 0; i < sz; ++i)
    {
        Vertex* otherEnd = GetOpposingSide(g, edges[i], v);
        depth_first_search(g, aaNID, otherEnd, id);
    }
}


void neuron_identification(Grid& g)
{
    // in case the IDs have already been calculated: do nothing
	static bool hasBeenDone = false;
	if (hasBeenDone) return;

    // synapse index attachment available?
    typedef Attachment<uint> ANeuronID;

    if (!GlobalAttachments::is_declared("neuronID"))
        UG_THROW("GlobalAttachment 'NeuronID' not declared.");

    ANeuronID aNID(GlobalAttachments::attachment<ANeuronID>("neuronID"));

    if (!g.has_vertex_attachment(aNID))
    	g.attach_to_vertices(aNID);

    Grid::VertexAttachmentAccessor<ANeuronID> aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(g, aNID);


    // initialize with -1
    VertexIterator it = g.begin<Vertex>();
    VertexIterator it_end = g.end<Vertex>();
    for (; it != it_end; ++it)
        aaNID[*it] = (uint) -1;

    uint nid = (uint) -1;
    for (it = g.begin<Vertex>(); it != it_end; ++it)
        if (aaNID[*it] == (uint) -1)
            depth_first_search(g, aaNID, *it, ++nid);

#ifdef UG_PARALLEL
    // in the parallel case, we need to offset distributed IDs
    unsigned long num_id = ++nid;
    if (pcl::NumProcs() > 1)
    {
    	size_t numProcs = (size_t) pcl::NumProcs();
    	int rank = pcl::ProcRank();
    	unsigned long* gatheredNIDcounts = NULL;
    	if (rank == 0)
    		gatheredNIDcounts = new unsigned long[numProcs];
    	pcl::ProcessCommunicator pc;
    	pc.gather(&num_id, 1, PCL_DT_UNSIGNED_LONG, gatheredNIDcounts, 1, PCL_DT_UNSIGNED_LONG, 0);

    	// calculate offsets
    	if (rank == 0)
    	{
    		unsigned long nextOffset = 0;
    		for (size_t i = 0; i < numProcs; ++i)
    		{
    			nextOffset += gatheredNIDcounts[i];
    			gatheredNIDcounts[i] = nextOffset - gatheredNIDcounts[i];
    		}
    	}

    	// scatter offsets
    	unsigned long offset;
    	pc.scatter(gatheredNIDcounts, 1, PCL_DT_UNSIGNED_LONG, &offset, 1, PCL_DT_UNSIGNED_LONG, 0);

		if (rank == 0)
			delete[] gatheredNIDcounts;

    	// add the offset to all ids
    	for (it = g.begin<Vertex>(); it != it_end; ++it)
			aaNID[*it] += (uint) offset;
    }
#endif

    // do not calculate IDs again
	hasBeenDone = true;
}



void save_neuron_to_swc
(
    std::string ugxFileName,
    size_t neuronIndex,
    std::string swcFileName,
    number scale
)
{
    Domain3d dom;
    try {LoadDomain(dom, ugxFileName.c_str());}
    UG_CATCH_THROW("Domain could not be loaded from file " << ugxFileName << ".");

    SmartPtr<MultiGrid> mg = dom.grid();
    Domain3d::position_accessor_type& aaPos = dom.position_accessor();
    Domain3d::subset_handler_type& sh = *dom.subset_handler();


    // get access to diameter attachment (or throw if not available)
    if (!GlobalAttachments::is_declared("diameter"))
        UG_THROW("GlobalAttachment 'diameter' not declared.");

    ANumber aDiam(GlobalAttachments::attachment<ANumber>("diameter"));
    //UG_COND_THROW(mg->has_vertex_attachment(aDiam), "Diameter attachment not available.");
    if (!mg->has_vertex_attachment(aDiam))
        mg->attach_to_vertices(aDiam);
    Grid::VertexAttachmentAccessor<ANumber> aaDiam(*mg, aDiam);


    // find a vertex of the neuron with required index
    typedef Attachment<uint> ANeuronID;

    if (!GlobalAttachments::is_declared("neuronID"))
        UG_THROW("GlobalAttachment 'NeuronID' not declared.");

    ANeuronID aNID(GlobalAttachments::attachment<ANeuronID>("neuronID"));
    Grid::VertexAttachmentAccessor<ANeuronID> m_aaNID;

    bool nidAvail = mg->has_vertex_attachment(aNID);

    Vertex* start = NULL;

    // if neuron indices are already available in the grid, use them
    if (nidAvail)
    {
        m_aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(*mg, aNID);

        // iterate vertices until correct nid is encountered
        VertexIterator it = mg->begin<Vertex>();
        VertexIterator it_end = mg->end<Vertex>();
        while (it != it_end && m_aaNID[*it] != neuronIndex)
            ++it;

        UG_COND_THROW(it == it_end, "A neuron with the required index " << neuronIndex
            << " is not contained in the grid " << ugxFileName << ".")

        start = *it;
    }

    // if neuron indices are not available, iterate through neurons until the required one is reached
    else
    {
        // iterate neurons
        mg->attach_to_vertices(aNID);

        Grid::VertexAttachmentAccessor<ANeuronID> aaNID(*mg, aNID);

        // initialize with -1
        VertexIterator it = mg->begin<Vertex>();
        VertexIterator it_end = mg->end<Vertex>();
        for (; it != it_end; ++it )
            aaNID[*it] = -1;

        uint nid = (uint) -1;
        for (it = mg->begin<Vertex>(); it != it_end; ++it)
        {
            if (aaNID[*it] == (uint) -1)
            {
                if (++nid == neuronIndex) break;
                depth_first_search(*mg, aaNID, *it, nid);
            }
        }

        UG_COND_THROW(nid < neuronIndex, "Neuron with index " << neuronIndex << " cannot be saved.\n"
            "Only " << nid+1 << " neurons contained in grid.");

        start = *it;
    }


    // analyze subset names to find out corresponding swc-types
    size_t nss = sh.num_subsets();
    std::vector<size_t> vType(nss);
    bool soma_subset_present = false;
    for (size_t i = 0; i < nss; ++i)
    {
        std::string name(sh.get_subset_name(i));
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name.find("SOMA") != std::string::npos)
        {
            soma_subset_present = true;
            vType[i] = 1;
        }
        else if (name.find("AXON") != std::string::npos)
            vType[i] = 2;
        else if (name.find("APIC") != std::string::npos)
            vType[i] = 4;
        else if (name.find("DEND") != std::string::npos)
            vType[i] = 3;
        else vType[i] = 0;
    }

    if (!soma_subset_present)
        UG_LOGN("Warning: No somatic subset could be identified.")

    // advance current vertex to soma (if identifiable)
    if (soma_subset_present)
    {
        mg->begin_marking();
        std::queue<Vertex*> q; // corresponds to breadth-first
        q.push(start);
        while (!q.empty())
        {
            Vertex* v = q.front();
            if (vType[sh.get_subset_index(v)] == 1) break;
            mg->mark(v);
            q.pop();

            // push neighboring elems to queue
            Grid::traits<Edge>::secure_container edges;
            mg->associated_elements(edges, v);

            size_t sz = edges.size();
            for (size_t e = 0; e < sz; ++e)
            {
                Vertex* otherEnd = GetOpposingSide(*mg, edges[e], v);
                if (!mg->is_marked(otherEnd))
                    q.push(otherEnd);
            }
        }
        mg->end_marking();

        if (q.empty())
            UG_LOGN("Warning: No soma vertex could be found in the requested neuron.")
        else
            start = q.front();
    }

    // write the neuron to file
    std::ofstream outFile(swcFileName.c_str(), std::ios::out);
    UG_COND_THROW(!outFile.is_open(), "Could not open output file '" << swcFileName << "'.");

    outFile << "# This file has been generated by UG4 from the original file "
            << ugxFileName << "." << std::endl;
    outFile << "# It contains neuron " << neuronIndex << " from that geometry." << std::endl;

    std::stack<std::pair<Vertex*, int> > stack; // corresponds to depth-first
    stack.push(std::make_pair(start, -1));

    mg->begin_marking();
    int ind = 0;   // by convention, swc starts with index 1
    bool all_types_identified = true;
    while (!stack.empty())
    {
        // get all infos regarding vertex
        std::pair<Vertex*, int>& info = stack.top();
        Vertex* v = info.first;
        int conn = info.second;
        stack.pop();

        // mark curr vrt
        mg->mark(v);

        size_t type = vType[sh.get_subset_index(v)];
        if (!type) all_types_identified = false;

        const Domain3d::position_type& coord = aaPos[v];

        number radius = 0.5*scale*aaDiam[v];

        // write line to file
        outFile << ++ind << " " << type << " "
            << coord[0]*scale << " " << coord[1]*scale << " " << coord[2]*scale << " "
            << radius << " " << conn << std::endl;

        // push neighboring elems to queue
        Grid::traits<Edge>::secure_container edges;
        mg->associated_elements(edges, v);

        size_t sz = edges.size();
        for (size_t e = 0; e < sz; ++e)
        {
            Vertex* otherEnd = GetOpposingSide(*mg, edges[e], v);
            if (!mg->is_marked(otherEnd))
                stack.push(std::make_pair(otherEnd, ind));
        }
    }
    mg->end_marking();

    if (!all_types_identified)
        UG_LOGN("WARNING: Some vertex type(s) - soma, dendrite, axon, etc. -\n"
            "could not be identified using the subset names.\n"
            << "To ensure correct types in the resulting swc file, the ugx subset names\n"
            "need to contain one of the strings \"SOMA\", \"AXON\", \"DEND\", \"APIC\",\n"
            "upper/lower case can be ignored.");

    outFile.close();
}




size_t innermost_neuron_id_in_subset(const std::string& ss, ConstSmartPtr<MGSubsetHandler> sh)
{
	// get subset id from name
	SubsetGroup ssg(sh, ss);
	UG_COND_THROW(ssg.size() != 1,
		"Could not determine a unique subset ID from given subset name '" << ss << "'.");
	int si = ssg[0];

	// position access
	const Domain3d::position_accessor_type& aaPos
		= Grid::VertexAttachmentAccessor<APosition>(*sh->grid(), aPosition);

	// neuron IDs
	typedef Attachment<uint> ANeuronID;

	if (!GlobalAttachments::is_declared("neuronID"))
		UG_THROW("GlobalAttachment 'NeuronID' not declared.");

	ANeuronID aNID(GlobalAttachments::attachment<ANeuronID>("neuronID"));
	Grid::VertexAttachmentAccessor<ANeuronID> m_aaNID;

	if (!sh->grid()->has_vertex_attachment(aNID))
	{
		sh->grid()->attach_to_vertices(aNID);
		neuron_identification(*sh->grid());
	}
	m_aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(*sh->grid(), aNID);


	// calculate center
	size_t n = 0;
	vector3 center(0.0);
	ConstVertexIterator it = sh->begin<Vertex>(si, 0);
	ConstVertexIterator it_end = sh->end<Vertex>(si, 0);
	for (; it != it_end; ++it)
	{
		++n;
		VecAdd(center, center, aaPos[*it]);
	}

	// communicate center and number
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		n = pc.allreduce(n, PCL_RO_SUM);
		vector3 globCenter;
		pc.allreduce(&center[0], &globCenter[0], 3, PCL_RO_SUM);
		center = globCenter;
	}
#endif

	VecScale(center, center, (number) 1.0 / n);

	// find vertex nearest to center
	uint minID = 0;
	number minSqDist = std::numeric_limits<number>::max();
	for (it = sh->begin<Vertex>(si, 0); it != it_end; ++it)
	{
		number distSq = VecDistanceSq(center, aaPos[*it]);
		if (distSq < minSqDist)
		{
			minSqDist = distSq;
			minID = m_aaNID[*it];
		}
	}

	// communicate
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		size_t nProc = pcl::NumProcs();

		number* minDists = NULL;
		if (pcl::ProcRank() == 0)
			minDists = new number[nProc];
		pc.gather(&minSqDist, 1, PCL_DT_DOUBLE, minDists, 1, PCL_DT_DOUBLE, 0);

		size_t minProc = 0;
		if (pcl::ProcRank() == 0)
		{
			for (size_t i = 1; i < nProc; ++i)
			{
				if (minDists[i] < minSqDist)
				{
					minSqDist = minDists[i];
					minProc = i;
				}
			}
			delete[] minDists;
		}
		pc.broadcast(minProc, 0);
		pc.broadcast(minID, minProc);
	}
#endif

	return minID;
}




number subset_length(int si, ConstSmartPtr<MGSubsetHandler> sh)
{
    // check domain, subset handler
    UG_COND_THROW(!sh.valid(), "Invalid subset handler.");

    const Domain3d::position_accessor_type& aaPos
        = Grid::VertexAttachmentAccessor<APosition>(*sh->grid(), aPosition);

    number totLength = 0.0;

    geometry_traits<Edge>::const_iterator it = sh->begin<Edge>(si, 0);
    geometry_traits<Edge>::const_iterator it_end = sh->end<Edge>(si, 0);
    for (; it != it_end; ++it)
        totLength += EdgeLength(*it, aaPos);

    return totLength;
}


number subset_length(const char* subset, ConstSmartPtr<MGSubsetHandler> sh)
{
    // check domain, subset handler
    UG_COND_THROW(!sh.valid(), "Invalid subset handler.");

    // get subset index from subset name
    int si = sh->get_subset_index(subset);

    // forward
    return subset_length(si, sh);
}



#ifdef UG_DIM_3
	template void scale_domain<Domain3d>(SmartPtr<Domain3d>, number);
	template bool is_acyclic<Domain3d>(SmartPtr<Domain3d>, int);
	template int check_presyn_indices<Domain3d>(SmartPtr<Domain3d>, int);
	template int check_domain<Domain3d>(SmartPtr<Domain3d>, int);
#endif


} // namespace cable_neruon
} // namespace ug



