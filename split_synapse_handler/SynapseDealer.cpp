/*
 * SynapseDealer.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#include "SynapseDealer.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

SynapseDealer::SynapseDealer() {
	// TODO tbd: Maybe register here?

}

SynapseDealer::~SynapseDealer() {
	//cleanup generator pointers
	for(std::map<std::string, ISynapseGenerator*>::iterator it = m_register.begin();
		it != m_register.end(); ++it) {
		delete it->second;
	}
}

void SynapseDealer::register_synapsetype(std::string t, ISynapseGenerator* image)
{
	m_register[t] = image;
}

IBaseSynapse* SynapseDealer::deal(std::string t)
{
	std::map<std::string, ISynapseGenerator*>::iterator it = m_register.find(t);
	//if t cannot be found return 0 else return the generated template of generate()
	return (it == m_register.end()) ? 0 : it->second->generate();
}

SynapseDealer *SynapseDealer::m_instance = 0;

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
