/*
 * order.cpp
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#include "order.h"

#include "common/error.h"
#include "lib_disc/domain.h"


namespace ug {
namespace cable {


/// help class to provide compare operator for indices based on their degree
/**
 * This class is used to provide an ordering for indices. The ordering relation
 * is based on the connectivity-degree, i.e. on the number of connections the
 * index has. The indices with less connections are ordered first.
 */
struct CompareDegree
{
	///	constructor, passing field with connections for each index
	CompareDegree
	(
		const std::vector<std::vector<size_t> >& vConn
	) : m_vCon(vConn){}

	///	comparison operator
	bool operator() (size_t i, size_t j)
	{
		UG_ASSERT(i < m_vCon.size() && j < m_vCon.size(), "size mismatch!");

		return (m_vCon[i].size() < m_vCon[j].size());
	}

private:
	///	storage field of connections of each index
	const std::vector<std::vector<size_t> >& m_vCon;
};

struct CompareIndex
{
	CompareIndex(const std::vector<size_t>& vInd)
	: m_vInd(vInd) {}

	bool operator() (size_t i, size_t j)
	{
		return m_vInd[i] < m_vInd[j];
	}

private:
	const std::vector<size_t>& m_vInd;
};


void compute_cuthillmckee_order
(
	std::vector<size_t>& vNewIndex,
	std::vector<std::vector<size_t> >& vvConnection
)
{
	PROFILE_FUNC();
	std::vector<size_t> vNewOrder;

	size_t conn_sz = vvConnection.size();

	// create flag list to remember already handled indices
	std::vector<bool> vHandled(conn_sz, false);

	// sort connections wrt. size of known neighbors
	std::vector<size_t> vConnSort(conn_sz);
	std::vector<size_t> vConnReverseSort(conn_sz);
	for (size_t c = 0; c < conn_sz; ++c)
		vConnSort[c] = vConnReverseSort[c] = c;

	CompareDegree cd(vvConnection);
	std::sort(vConnSort.begin(), vConnSort.end(), cd);
	CompareIndex ci(vConnSort);
	std::sort(vConnReverseSort.begin(), vConnReverseSort.end(), ci);

	// vConnSort now contains the indices of the vvConnection vector
	// sorted by the number of entries
	// vConnReverseSort contains the index an entry of vvConnection has
	// in vConnSort (inverse permutation)

	// sort adjacent index by degree
	for (size_t i = 0; i < conn_sz; ++i)
	{
		// indices with no adjacent indices are marked as handled (and skipped)
		if (vvConnection[vConnSort[i]].size() == 0)
			vHandled[i] = true;
		else
			std::sort(vvConnection[i].begin(), vvConnection[i].end(), cd);
	}

	size_t start = 0;
	while(true)
	{
		// find first unhandled index
		size_t i_notHandled = start;
		for (; i_notHandled < vHandled.size(); ++i_notHandled)
			if (!vHandled[i_notHandled]) {start = i_notHandled; break;}

		// check if any unhandled vertex left
		if (i_notHandled == vHandled.size()) break;

		// add start vertex to mapping
		vNewOrder.push_back(vConnSort[start]);
		vHandled[start] = true;

		// create queue of adjacent vertices
		std::queue<size_t> qAdjacent;
		for (size_t i = 0; i < vvConnection[vConnSort[start]].size(); ++i)
		{
			const size_t ind = vConnReverseSort[vvConnection[vConnSort[start]][i]];

			if (!vHandled[ind] && ind != start)
				qAdjacent.push(ind);
		}

		// work off queue
		while (!qAdjacent.empty())
		{
			const size_t front = qAdjacent.front();

			if (!vHandled[front])
			{
				// add to mapping
				vNewOrder.push_back(vConnSort[front]);
				vHandled[front] = true;

				// add adjacent vertices to queue
				for (size_t i = 0; i < vvConnection[vConnSort[front]].size(); ++i)
				{
					const size_t ind = vConnReverseSort[vvConnection[vConnSort[front]][i]];

					if (!vHandled[ind] && ind != front)
						qAdjacent.push(ind);
				}
			}

			qAdjacent.pop();
		}
	}

	// create list of mapping
	vNewIndex.clear();
	vNewIndex.resize(conn_sz, (size_t)-1);

	// write new indices into out array
	size_t cnt = 0;
	for (size_t oldInd = 0; oldInd < conn_sz; ++oldInd)
	{
		// skip non-sorted indices
		if (vvConnection[oldInd].size() == 0) continue;

		// get old index
		UG_ASSERT(cnt < vNewOrder.size(), "cnt: "<<cnt<<", ordered: "<<vNewOrder.size())
		const size_t newInd = vNewOrder[vNewOrder.size() - 1 - cnt]; ++cnt;
		UG_ASSERT(newInd < vNewIndex.size(), "newInd: "<<newInd<<", size: "<<vNewIndex.size())

		// set new index to order
		vNewIndex[newInd] = oldInd;
	}

	// check if all ordered indices have been written
	if (cnt != vNewOrder.size())
		UG_THROW("OrderCuthillMcKee: Not all indices sorted that must be sorted: "
				<< cnt << " written, but should write: " << vNewOrder.size());

//	fill non-sorted indices
/*	TODO: This is definitely wrong in general if DoFs are allowed not to have any connections!
 *	Suppose, we have N nodes containing one DoF each (indexed 0, 1, ... N-1).
 *	Now, DoF i might not have any connections, while all the other DoFs do. In that case,
 *	the only possible value for vNewIndex[i] is i, as all the other values are already taken
 *	by the other indices. But in general, vNewIndex[i-1] will not be i-1 as it would have to
 *	be in order to ensure vNewIndex[i] = i!
 *	Furthermore, with bad luck, DoF 0 might be disconnected. In that case, this loop does not
 *	even set _any_ value and vNewIndex[0] will contain -1.
 *
 *	BEWARE, however: One cannot simply assign vNewIndex[i] = i (as would seem to be a logical
 *	choice of ordering), since DoF i-1 and DoF i might belong to the same node (in a case where
 *	there is more than one unknown per node) and must therefore have consecutive indices!
 *
 *	The intended usage of this method is like this:
 *	vvConnection is filled by calling dofDistr.get_connections(vvConnection).
 *	Then for any geometry that does not contain isolated (i.e. unconnected) vertices,
 *	it will contain connections EXACTLY for every first DoFs of any node i.e.:
 *	- It does not contain any connections for any other DoF BUT the first of every node.
 *	- There is not any node for which the first DoF has no connections.
 *
 *	If these conditions are satisfied the code will work.
 *	Otherwise (like in the above example) this will not be the case in general.
 *
 *	Best solution at the moment: Avoid DoFs that are disconnected!
 */
	for (size_t i = 1; i < vNewIndex.size(); ++i)
		if (vNewIndex[i] == (size_t)-1)
			vNewIndex[i] = vNewIndex[i-1] + 1;

	//CheckPermutationBijective(vNewIndex);
}



void order_cuthillmckee(DoFDistribution& dofDistr)
{
	// get adjacency graph
	std::vector<std::vector<size_t> > vvConnection;
	try
	{
		dofDistr.get_connections(vvConnection);
	}
	UG_CATCH_THROW("OrderCuthillMcKee: No adjacency graph available.");

	// get mapping for cuthill-mckee order
	std::vector<size_t> vNewIndex;
	compute_cuthillmckee_order(vNewIndex, vvConnection);

	//reorder indices
	dofDistr.permute_indices(vNewIndex);
}



template <typename TDomain>
void order_cuthillmckee(ApproximationSpace<TDomain>& approxSpace)
{
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	for(size_t i = 0; i < vDD.size(); ++i)
		order_cuthillmckee(*vDD[i]);
}



#ifdef UG_DIM_1
	template void order_cuthillmckee<Domain1d>(ApproximationSpace<Domain1d>& approxSpace);
#endif
#ifdef UG_DIM_2
	template void order_cuthillmckee<Domain2d>(ApproximationSpace<Domain2d>& approxSpace);
#endif
#ifdef UG_DIM_3
	template void order_cuthillmckee<Domain3d>(ApproximationSpace<Domain3d>& approxSpace);
#endif


} // namespace cable
} // namespace ug


