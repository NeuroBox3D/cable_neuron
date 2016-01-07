/*
 * synapse_info_attachment_handler.h
 *
 *  Created on: 20.05.2015
 *      Author: mbreit
 */
#ifndef PLUGINS__SYNAPSE_HANDLER__SYNAPSE_INFO_ATTACHMENT_HANDLER_H_
#define PLUGINS__SYNAPSE_HANDLER__SYNAPSE_INFO_ATTACHMENT_HANDLER_H_


#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"
#include "../grid/synapse_info.h"		// synapse_traits

#include <vector>

namespace ug{
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
} // end namespace ug

#endif // PLUGINS__SYNAPSE_HANDLER__SYNAPSE_INFO_ATTACHMENT_HANDLER_H_
