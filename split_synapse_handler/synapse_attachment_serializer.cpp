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
	size_t sz = m_aa[o].size();
	Serialize(out, sz);
	for (size_t i = 0; i < sz; ++i)
	{
		// save synapse type to buffer
		Serialize(out, m_aa[o][i]->name());

		// save synapse to buffer
		Serialize(out, *m_aa[o][i]);

		// delete object
		delete m_aa[o][i];
	}
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
			std::string name = Deserialize<std::string>(in);

			// create new synapse
			IBaseSynapse* s = SynapseDealer::instance()->deal(name);
			UG_LOGN("  synapse type: " << name);

			// read values from buffer
			Deserialize(in, s);

			// save pointer in attachment
			m_aa[o][i] = s;
		}
	}
	UG_CATCH_THROW("Synapse attachment distribution did not succeed.")
}


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
