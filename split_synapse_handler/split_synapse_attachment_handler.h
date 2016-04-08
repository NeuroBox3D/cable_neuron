/*
 * splitsynapseattachmenthandler.h
 *
 *  Created on: Apr 6, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_ATTACHMENT_HANDLER_H_
#define SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_ATTACHMENT_HANDLER_H_

#include <vector>
#include "base_synapse.h"

#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/copy_attachment_handler.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class SplitSynapseAttachmentHandler
: public CopyAttachmentHandler<Edge, Attachment<std::vector<IBaseSynapse*> > >
{
public:
	SplitSynapseAttachmentHandler() {}
	virtual ~SplitSynapseAttachmentHandler() {}
protected:
	virtual void copy(Edge* parent, Edge* child);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SPLIT_SYNAPSE_ATTACHMENT_HANDLER_H_ */
