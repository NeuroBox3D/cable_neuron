/*
 * synapseattachmenthandler.cpp
 *
 *  Created on: Apr 6, 2016
 *      Author: lreinhardt
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
