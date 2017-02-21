/*
 * synapse_info_attachment_handler.h
 *
 *  Created on: 20.05.2015
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__GRID__SYNAPSE_INFO_ATTACHMENT_HANDLER_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__GRID__SYNAPSE_INFO_ATTACHMENT_HANDLER_H__


#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"
#include "../grid/synapse_info.h"		// synapse_traits

#include <vector>

namespace ug{
namespace cable_neuron {
namespace synapse_handler{


/**
 * @brief handler for SynapseInfo attachments in a multi-grid
 *
 * This class implements an attachment handler for the synapse attachment used
 * in the NETISynapseHandle class.
 *
 * Synapse edge attachments need to be propagated not merely by copying as synapses
 * are point processes on edges and thus only "rightfully" belong to one of the child
 * edges in the multi-grid hierarchy.
 *
 * This attachment handler derives from CopyAttachmentHandler and overloads the
 * CopyAttachmentHandler<Edge, Attachment<std::vector<SynapseInfo> > >::copy()
 * method in order to be able to copy only to the correct child edge.
 *
 */
class SynapseInfoAttachmentHandler
: public CopyAttachmentHandler<Edge, Attachment<std::vector<SynapseInfo> > >
{
	public:
		typedef synapse_traits<> STV;

		/// constructor
		SynapseInfoAttachmentHandler() {};

		/// destructor
		virtual ~SynapseInfoAttachmentHandler() {};

	protected:
		virtual void copy(Edge* parent, Edge* child);
};

} // end namespace synapse_handler
} // namespace cable_neuron
} // end namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__GRID__SYNAPSE_INFO_ATTACHMENT_HANDLER_H__
