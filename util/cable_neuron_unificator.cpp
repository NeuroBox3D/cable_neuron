/*
 * cable_neuron_unificator.cpp
 *
 *  Created on: 2019-03-22
 *      Author: mbreit
 */

#include "cable_neuron_unificator.h"


namespace ug {
namespace cable_neuron {


CableNeuronUnificator::CableNeuronUnificator()
{}


CableNeuronUnificator::~CableNeuronUnificator()
{}


void CableNeuronUnificator::unify
(
	MultiGrid* mg,
	int lvl,
	int localOffset,
	const Grid::AttachmentAccessor<Edge, AElemIndex>& aaElemInd, // local indices!
	const Grid::AttachmentAccessor<side_t, AElemIndices>& aaSideElemInd, // global indices!
	std::vector<std::pair<int, int> >& unificationPairs // global indices!
) const
{
	Grid::traits<Edge>::secure_container el;

	typedef geometry_traits<Edge>::const_iterator const_edge_iter;
	const_edge_iter it = mg->begin<Edge>(lvl);
	const_edge_iter itEnd = mg->end<Edge>(lvl);
	for (; it != itEnd; ++it)
	{
		Edge* e = *it;
		int ind = aaElemInd[e];

		// ignore ghosts
		if (ind == -1)
			continue;

		for (size_t i = 0; i < 2; ++i)
		{
			mg->associated_elements(el, e->vertex(i));
			const size_t elSz = el.size();
			for (size_t e2 = 0; e2 < elSz; ++e2)
			{
				Edge* other = el[e2];

				if (other <= e)
					continue;

				int ind2 = aaElemInd[other];
				UG_COND_THROW(ind2 == -1, "CableNeuronUnificator is supposed to hold neurons together,\n"
					"but this neuron is already cut by the partitioning.");

				unificationPairs.push_back(std::make_pair(ind + localOffset, ind2 + localOffset));
			}
		}
	}
}


}  // namespace cable_neuron
}  // namespace ug
