/*
 * SynapseDealer.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "synapse_dealer.h"
#include "common/log.h"
#include "common/assert.h"
#include "pcl/pcl_process_communicator.h"


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
	std::map<std::string, size_t>::const_iterator it = m_register.find(name);
	UG_COND_THROW(it == m_register.end(), "Requested unique ID for synapse type named '"
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
		recBufOr = new char[sz];
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

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
