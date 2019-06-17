/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-03-22
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
