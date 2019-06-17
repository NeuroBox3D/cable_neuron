/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Lukas Reinhardt
 * Creation date: 2016-04-06
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

#include "synapse_attachment_handler.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

void SynapseAttachmentHandler::copy(Edge* parent, Edge* child)
{
	// remove anything from child attachment
	m_aa[child].resize(0);

	// iterate over all synapses of the parent edge
	for (size_t j = 0; j < m_aa[parent].size(); ++j)
	{
		IBaseSynapse* s = m_aa[parent][j];

		// find out whether child is parent's left or right child
		Vertex* leftChildVrt = (*child)[0];
		Vertex* leftParentVrt = (*parent)[0];

		// left child
		if (m_spMG->get_child_vertex(leftParentVrt) == leftChildVrt)
		{
			// synapse belongs to other child
			if (s->location() >= 0.5)
				continue;
			// synapse belongs to child
			else
			{
				// copy info
				m_aa[child].push_back(s);

				// adapt local coord
				m_aa[child].back()->set_location(2.0*s->location());
			}
		}
		// right child
		else if (m_spMG->get_child_vertex((*parent)[1]) == (*child)[1])
		{
			// synapse belongs to other child
			if (s->location() < 0.5)
				continue;
			// synapse belongs to child
			else
			{
				// copy info
				m_aa[child].push_back(s);

				// adapt local coord
				m_aa[child].back()->set_location(2.0*s->location() - 1.0);
			}
		}
		// error
		else
			UG_THROW("Child edge is neither left nor right child.");
	}
}

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
