/*
 * SynapseDealer.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "synapse_dealer.h"
#include "common/log.h"

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


size_t SynapseDealer::unique_id(std::string name)
{
	std::map<std::string, size_t>::const_iterator it = m_register.find(name);
	UG_COND_THROW(it == m_register.end(), "Requested unique ID for synapse type named '"
		<< name << "' which is not registered.");

	return it->second;
}


IBaseSynapse* SynapseDealer::deal(std::string name)
{
	// find index for synapse type
	size_t uid;
	try {uid = unique_id(name);}
	UG_CATCH_THROW("Requested synapse type '" << name << "' which is not registered.");

	// get generator using index and generate synapse
	return m_vSynGen[uid]->generate();
}

IBaseSynapse* SynapseDealer::deal(size_t unique_id)
{
	return m_vSynGen[unique_id]->generate();
}

size_t SynapseDealer::size_of(std::string name)
{
	// find index for synapse type
	size_t uid;
	try {uid = unique_id(name);}
	UG_CATCH_THROW("Requested synapse type '" << name << "' which is not registered.");

	return m_vSize[uid];
}

size_t SynapseDealer::size_of(size_t unique_id)
{
	return m_vSize[unique_id];
}

size_t SynapseDealer::get_unique_id()
{
	static size_t uid = 0;
	return uid++;
}



SynapseDealer *SynapseDealer::m_instance = 0;

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
