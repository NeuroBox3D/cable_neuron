/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Lukas Reinhardt
 * Creation date: 2016-03-22
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

#include "synapse_dealer.h"
#include "common/log.h"
#include "common/assert.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

SynapseDealer::SynapseDealer()
{}

SynapseDealer::~SynapseDealer()
{
	// cleanup generator pointers
	size_t sz = m_vSynGen.size();
	for (size_t i = 0; i < sz; ++i)
		delete m_vSynGen[i];
}


size_t SynapseDealer::unique_id(const std::string& name)
{
	std::map<std::string, size_t>::const_iterator it = m_synReg.find(name);
	UG_COND_THROW(it == m_synReg.end(), "Requested unique ID for synapse type named '"
		<< name << "' which is not registered.");

	return it->second;
}


IBaseSynapse* SynapseDealer::deal(const std::string& name)
{
	// find index for synapse type
	size_t uid;
	try {uid = unique_id(name);}
	UG_CATCH_THROW("Requested synapse type '" << name << "' which is not registered.");

	// get generator using index and generate synapse
	return deal(uid);
}

IBaseSynapse* SynapseDealer::deal(size_t unique_id)
{
	UG_ASSERT(m_vSynGen.size() > unique_id,
			  "Requested synapse of unknown ID " << unique_id << ".");
	return m_vSynGen[unique_id]->generate();
}

size_t SynapseDealer::size_of(const std::string& name)
{
	// find index for synapse type
	size_t uid;
	try {uid = unique_id(name);}
	UG_CATCH_THROW("Requested synapse type '" << name << "' which is not registered.");

	return size_of(uid);
}

size_t SynapseDealer::size_of(size_t unique_id)
{
	UG_ASSERT(m_vSynGen.size() > unique_id,
				  "Requested synapse type size of unknown ID " << unique_id << ".");
	return m_vSize[unique_id];
}


const std::vector<size_t>& SynapseDealer::non_data_bytes(size_t uid) const
{
	UG_ASSERT(m_vNoOverwriteBytes.size() > uid,
				  "Requested non-data bytes for synapse of unknown ID " << uid << ".");
	return m_vNoOverwriteBytes[uid];
}



size_t SynapseDealer::get_unique_id()
{
	static size_t uid = 0;
	return uid++;
}


void SynapseDealer::find_no_overwrite_bytes
(
	const char* buf,
	size_t uid
)
{
#ifdef UG_PARALLEL
	if (pcl::NumProcs() < 2) return;

	size_t sz = m_vSize[uid];

	pcl::ProcessCommunicator com;

	char* recBufAnd = NULL;
	char* recBufOr = NULL;
	if (pcl::ProcRank() == 0)
	{
		recBufAnd = new char[sz];
		std::memset(recBufAnd, 0, sz);
		recBufOr = new char[sz];
		std::memset(recBufOr, 0, sz);
	}

	com.reduce(buf, recBufAnd, sz, PCL_DT_CHAR, PCL_RO_BAND, 0);
	com.reduce(buf, recBufOr, sz, PCL_DT_CHAR, PCL_RO_BOR, 0);

	std::vector<size_t>& vNOB = m_vNoOverwriteBytes[uid];

	if (pcl::ProcRank() == 0)
	{
		for (size_t b = 0; b < sz; ++b)
		{
			// if not all bits are equal on all procs, this is a no-overwrite byte
			if (~(recBufAnd[b] | ~recBufOr[b]))
				vNOB.push_back(b);
		}
		delete[] recBufAnd;
		delete[] recBufOr;
	}

	size_t nobSz = vNOB.size();

	com.broadcast(nobSz, 0);
	vNOB.resize(nobSz);
	if (nobSz)
		com.broadcast(&vNOB[0], nobSz, 0);

	/* DEBUG
	UG_LOG(uid << ":  ");
	for (size_t i = 0; i < nobSz; ++i)
		UG_LOG(vNOB[i] << " ")
	UG_LOGN("");
	*/
#endif
}


SynapseDealer *SynapseDealer::m_instance = 0;

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
