/*
 * synapse_info_attachment_handler.cpp
 *
 *  Created on: 20.05.2015
 *      Author: mbreit
 */


#include "synapse_info_attachment_handler.h"

namespace ug{
namespace synapse_handler{


void SynapseInfoAttachmentHandler::copy(Edge* parent, Edge* child)
{
	// remove anything from child attachment
	m_aa[child].resize(0);

	// iterate over all synapses of the parent edge
	for (size_t j = 0; j < m_aa[parent].size(); ++j)
	{
		SynapseInfo& info = m_aa[parent][j];

		// find out whether child is parent's left or right child
		Vertex* leftChildVrt = (*child)[0];
		Vertex* leftParentVrt = (*parent)[0];

		// left child
		if (m_spMG->get_child_vertex(leftParentVrt) == leftChildVrt)
		{
			// synapse belongs to other child
			if (STV::loc_coord(info) >= 0.5)
				continue;
			// synapse belongs to child
			else
			{
				// copy info
				m_aa[child].push_back(info);

				// adapt local coord
				STV::loc_coord(m_aa[child].back()) = 2.0*STV::loc_coord(info);
			}
		}
		// right child
		else if (m_spMG->get_child_vertex((*parent)[1]) == (*child)[1])
		{
			// synapse belongs to other child
			if (STV::loc_coord(info) < 0.5)
				continue;
			// synapse belongs to child
			else
			{
				// copy info
				m_aa[child].push_back(info);

				// adapt local coord
				STV::loc_coord(m_aa[child].back()) = 2.0*STV::loc_coord(info) - 1.0;
			}
		}
		// error
		else
			UG_THROW("Child edge is neither left nor right child.");
	}
}

} // end namespace synapse_handler
} // end namespace ug

