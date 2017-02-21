/*
 * synapse_attachment_handler.h
 *
 *  Created on: Apr 6, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_HANDLER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_HANDLER_H

#include "synapses/base_synapse.h"
#include <vector>

#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class SynapseAttachmentHandler
: public CopyAttachmentHandler<Edge, Attachment<std::vector<IBaseSynapse*> > >
{
public:
	SynapseAttachmentHandler() {}
	virtual ~SynapseAttachmentHandler() {}
protected:
	virtual void copy(Edge* parent, Edge* child);
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_HANDLER_H
