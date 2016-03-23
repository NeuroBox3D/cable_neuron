/*
 * SynapseDealer.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSEDEALER_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSEDEALER_H_

#include "IBaseSynapse.h"
#include <map>
#include <string>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class ISynapseGenerator
{
public:
	virtual IBaseSynapse* generate() = 0;
	virtual ~ISynapseGenerator() {}
};

template<typename TSyn>
class SynapseGenerator : public ISynapseGenerator
{
public:
	IBaseSynapse* generate() {return new TSyn;}
};

class SynapseDealer {
private:
	std::map<std::string, ISynapseGenerator*> m_register;
	static SynapseDealer *m_instance;
public:
	SynapseDealer();
	virtual ~SynapseDealer();

	void register_synapsetype(std::string t, ISynapseGenerator* image);
	IBaseSynapse* deal(std::string t);

	static SynapseDealer* instance() {
		if(m_instance == 0) {
			m_instance = new SynapseDealer;
		}
		return m_instance;
	}
};

SynapseDealer *SynapseDealer::m_instance = 0;

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSEDEALER_H_ */
