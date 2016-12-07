/*
 * synapse_attachment_serializer.cpp
 *
 *  Created on: 02.12.2016
 *      Author: mbreit
 */

#include "synapse_attachment_serializer.h"
#include "synapse_dealer.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

SynapseAttachmentSerializer::SynapseAttachmentSerializer(Grid& g, SynAttach& a)
: m_aa(g, a, true) {};

SynapseAttachmentSerializer::~SynapseAttachmentSerializer()
{}


void SynapseAttachmentSerializer::write_data(BinaryBuffer& out, Edge* o) const
{
	try
	{
		size_t sz = m_aa[o].size();
		Serialize(out, sz);
		for (size_t i = 0; i < sz; ++i)
		{
			// get unique ID for synapse type
			std::string name = m_aa[o][i]->name();
			SynapseDealer* synReg = SynapseDealer::instance();
			size_t uid = synReg->unique_id(name);

			// save synapse type to buffer
			Serialize(out, uid);

			// save synapse to buffer
			size_t sz = synReg->size_of(uid);
			out.write(reinterpret_cast<char*>(m_aa[o][i]), sz);

			// delete object
			delete m_aa[o][i];
		}
	}
	UG_CATCH_THROW("Something went wrong during serialization of synapse attachment.");
}


void SynapseAttachmentSerializer::read_data(BinaryBuffer& in, Edge* o)
{
	try
	{
		// read number of synapses
		size_t sz = Deserialize<size_t>(in);
		UG_LOGN("vector size: " << sz);
		m_aa[o].resize(sz);

		for (size_t i = 0; i < sz; ++i)
		{
			// read synapse type
			size_t uid = Deserialize<size_t>(in);

			// get synapse registry/dealer
			SynapseDealer* synReg = SynapseDealer::instance();

			// create new synapse
			IBaseSynapse* s = synReg->deal(uid);
			UG_LOGN("  synapse type: " << uid);

			// read values from buffer
			in.read(reinterpret_cast<char*>(s), synReg->size_of(uid));

			// save pointer in attachment
			m_aa[o][i] = s;
		}
	}
	UG_CATCH_THROW("Synapse attachment distribution did not succeed.")
}


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
