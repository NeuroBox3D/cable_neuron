/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-02
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

			// save synapse to buffer
			size_t sz = synReg->size_of(uid);
			size_t nonDataPosPos = 0;
			const std::vector<size_t>& nonDataBytes = synReg->non_data_bytes(uid);
			size_t nonDataPos = nonDataBytes.size() ? nonDataBytes[0] : sz;
			for (size_t j = 0; j < sz; ++j)
			{
				if (j != nonDataPos)
					out.write(reinterpret_cast<char*>(m_aa[e][i])+j, 1);
				else
					nonDataPos = nonDataBytes.size() > ++nonDataPosPos ? nonDataBytes[nonDataPosPos] : sz;
			}

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

			// read values from buffer
			size_t sz = synReg->size_of(uid);
			size_t nonDataPosPos = 0;
			const std::vector<size_t>& nonDataBytes = synReg->non_data_bytes(uid);
			size_t nonDataPos = nonDataBytes.size() ? nonDataBytes[0] : sz;
			for (size_t j = 0; j < sz; ++j)
			{
				if (j != nonDataPos)
					in.read(reinterpret_cast<char*>(s)+j, 1);
				else
					nonDataPos = nonDataBytes.size() > ++nonDataPosPos ? nonDataBytes[nonDataPosPos] : sz;
			}

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
