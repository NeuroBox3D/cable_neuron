/*
 * SynapseDealer_impl.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "common/error.h"
#include "common/log.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {


template <typename TSyn>
void SynapseDealer::register_synapse_type(const std::string& t)
{
	// check if name already exists and throw if that is the case
	std::map<std::string, size_t>::const_iterator it = m_register.find(t);
	UG_COND_THROW(it != m_register.end(), "Tried to register synapse type named "
			      "'" << t << "', but this name already exists.");

	// assign a unique id to each synapse type;
	// save corresponding generators and object sizes
	size_t uid = get_unique_id();
	ISynapseGenerator* sg = new SynapseGenerator<TSyn>();
	m_register[t] = uid;
	m_vSynGen.resize(uid+1);
	m_vSynGen[uid] = sg;
	m_vSize.resize(uid+1);
	m_vSize[uid] = sizeof(TSyn);

	// compute and save byte positions not to be overwritten
	// when communicating this object type from one proc to another
	m_vNoOverwriteBytes.resize(uid+1);

#ifdef UG_PARALLEL
	char* buf = new char[sizeof(TSyn)];
	memset(buf, 0, sizeof(TSyn));
	new (buf) TSyn();
	find_no_overwrite_bytes(buf, uid);
	delete[] buf;
#endif

	/* DEBUGGING
	UG_LOGN("Current synapse registry");
	for (std::map<std::string, size_t>::const_iterator it = m_register.begin(); it != m_register.end(); ++it)
		UG_LOGN(it->first << ": " << it->second);
	UG_LOGN("");
	*/
}


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
