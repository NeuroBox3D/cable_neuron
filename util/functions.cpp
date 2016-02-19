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
#include "pcl/pcl_base.h"


namespace ug {
namespace cable_neuron {


template <typename TDomain>
bool is_acyclic(SmartPtr<TDomain> dom)
{
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


#ifdef UG_DIM_1
	template bool is_acyclic<Domain1d>(SmartPtr<Domain1d>);
#endif
#ifdef UG_DIM_2
	template bool is_acyclic<Domain2d>(SmartPtr<Domain2d>);
#endif
#ifdef UG_DIM_3
	template bool is_acyclic<Domain3d>(SmartPtr<Domain3d>);
#endif


} // namespace cable_neruon
} // namespace ug



