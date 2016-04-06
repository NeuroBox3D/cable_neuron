/*
 * SynapseDealer.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_

#include <map>
#include <string>

#include "base_synapse.h"

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
	IBaseSynapse* generate() {return new TSyn();}
};

class SynapseDealer {
private:
	std::map<std::string, ISynapseGenerator*> m_register;
	static SynapseDealer *m_instance;
	~SynapseDealer();
	SynapseDealer();
public:

	template <typename TSyn> void register_synapse_type(std::string t);
	IBaseSynapse* deal(std::string t);

	static SynapseDealer* instance() {
		if(m_instance == 0) {
			m_instance = new SynapseDealer;
		}
		return m_instance;
	}
};


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */


#include "synapse_dealer_impl.h"

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSE_DEALER_H_ */
