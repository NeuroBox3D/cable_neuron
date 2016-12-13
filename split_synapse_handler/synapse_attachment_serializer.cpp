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
: m_aa(g, a, true), m_blubb(false) {};

SynapseAttachmentSerializer::~SynapseAttachmentSerializer()
{}


void SynapseAttachmentSerializer::write_data(BinaryBuffer& out, Edge* e) const
{
	try
	{
		SynapseDealer* synReg = SynapseDealer::instance();

		size_t sz = m_aa[e].size();
		Serialize(out, sz);
		for (size_t i = 0; i < sz; ++i)
		{
			// get unique ID for synapse type
			const std::string& name = m_aa[e][i]->name();
			size_t uid = synReg->unique_id(name);

			// save synapse type to buffer
			Serialize(out, uid);
/*
			// save synapse to buffer
			size_t sz = synReg->size_of(uid);
			size_t nonDataPosPos = 0;
			const std::vector<size_t>& nonDataBytes = synReg->non_data_bytes(uid);
			size_t nonDataPos = nonDataBytes.size() ? nonDataBytes[0] : sz;
			for (size_t j = 0; j < sz; ++j)
			{
				if (j != nonDataPos)
					;//out.write(reinterpret_cast<char*>(m_aa[o][i])+j, 1);
				else
					nonDataPos = nonDataBytes.size() > ++nonDataPosPos ? nonDataBytes[nonDataPosPos] : sz;
			}
*/
			// delete object
			delete m_aa[e][i];
		}
	}
	UG_CATCH_THROW("Something went wrong during serialization of synapse attachment.");
}


void SynapseAttachmentSerializer::read_data(BinaryBuffer& in, Edge* e)
{
	try
	{
		// read number of synapses
		size_t sz = Deserialize<size_t>(in);
		m_aa[e].resize(sz);

		// get synapse registry/dealer
		SynapseDealer* synReg = SynapseDealer::instance();

		for (size_t i = 0; i < sz; ++i)
		{
			// read synapse type
			size_t uid = Deserialize<size_t>(in);

			// create new synapse
			IBaseSynapse* s = synReg->deal(uid);
			//UG_LOGN("  synapse type: " << uid << "   (" << (size_t*) s << ")");

			/*
			UG_LOGN("synapse before reading from buffer: ")
			for (size_t j = 0; j < synReg->size_of(uid); ++j)
			{
				char* b = ((char*) s) + j;
				UG_LOG(std::hex << std::setfill('0') << std::setw(2) << (int) (*b & 0xFF) << " ");
			}
			UG_LOGN("");
			*/
/*
			// read values from buffer
			size_t sz = synReg->size_of(uid);
			size_t nonDataPosPos = 0;
			const std::vector<size_t>& nonDataBytes = synReg->non_data_bytes(uid);
			size_t nonDataPos = nonDataBytes.size() ? nonDataBytes[0] : sz;
			for (size_t j = 0; j < sz; ++j)
			{
				if (j != nonDataPos)
					;//in.read(reinterpret_cast<char*>(s)+j, 1);
				else
					nonDataPos = nonDataBytes.size() > ++nonDataPosPos ? nonDataBytes[nonDataPosPos] : sz;
			}
*/
			/*
			UG_LOGN("synapse after reading from buffer:  ")
			for (size_t j = 0; j < synReg->size_of(uid); ++j)
			{
				char* b = ((char*) s) + j;
				UG_LOG(std::hex << std::setfill('0') << std::setw(2) << (int) (*b & 0xFF) << " ");
			}
			UG_LOGN("");
			*/

			// save pointer in attachment
			m_aa[e][i] = s;
		}
	}
	UG_CATCH_THROW("Synapse attachment distribution did not succeed.");
}


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
