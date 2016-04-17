/*
 * functions.cpp
 *
 *  Created on: 17.02.2016
 *      Author: mbreit
 */


#include "functions.h"

#include <stack>
#include <utility>

#include "lib_disc/domain.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_grid/grid/grid_base_object_traits.h"
#include "lib_grid/global_attachments.h" // global attachments
#include "lib_grid/tools/surface_view.h"
#include "pcl/pcl_base.h"


namespace ug {
namespace cable_neuron {



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




#ifdef UG_DIM_1
	template bool is_acyclic<Domain1d>(SmartPtr<Domain1d>, int);
	template int check_presyn_indices<Domain1d>(SmartPtr<Domain1d>, int);
	template int check_domain<Domain1d>(SmartPtr<Domain1d>, int);
#endif
#ifdef UG_DIM_2
	template bool is_acyclic<Domain2d>(SmartPtr<Domain2d>, int);
	template int check_presyn_indices<Domain2d>(SmartPtr<Domain2d>, int);
	template int check_domain<Domain2d>(SmartPtr<Domain2d>, int);
#endif
#ifdef UG_DIM_3
	template bool is_acyclic<Domain3d>(SmartPtr<Domain3d>, int);
	template int check_presyn_indices<Domain3d>(SmartPtr<Domain3d>, int);
	template int check_domain<Domain3d>(SmartPtr<Domain3d>, int);
#endif


} // namespace cable_neruon
} // namespace ug



