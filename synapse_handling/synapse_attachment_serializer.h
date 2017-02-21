/*
 * synapse_attachment_serializer.h
 *
 *  Created on: 02.12.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H

#include "lib_grid/algorithms/serialization.h"   // GeomObjDataSerializer
#include "synapses/base_synapse.h"  // IBaseSynapse

#include <vector>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


/**
 * As the synapse attachment is a vector only containing pointers
 * to the actual synapse objects, we need to manually take care of communicating
 * these objects during domain decomposition and redistribution.
 */
class SynapseAttachmentSerializer
: public GeomObjDataSerializer<Edge>
{
	public:
		typedef Attachment<std::vector<IBaseSynapse*> > SynAttach;

	public:
		SynapseAttachmentSerializer(Grid& g, SynAttach& a);

		virtual ~SynapseAttachmentSerializer();

		/**
		 * @brief Write synapses of a specific edge to binary buffer.
		 * When synapse is written to buffer it is guaranteed to be moved to another process.
		 * This is why we need to delete it here via its pointer on the current proc.
		 */
		virtual void write_data(BinaryBuffer& out, Edge* o) const;

		/**
		 * @brief Read synapses of a specific edge from binary buffer.
		 * The data received in the buffer will trigger the creation of a new synapse
		 * and storing of a pointer to it in the attachment.
		 */
		virtual void read_data(BinaryBuffer& in, Edge* o);

	protected:
		Grid::AttachmentAccessor<Edge, SynAttach> m_aa;

		bool m_blubb;
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_ATTACHMENT_SERIALIZER_H
