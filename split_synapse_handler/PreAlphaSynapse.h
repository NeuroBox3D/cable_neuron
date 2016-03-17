/*
 * PreAlphaSynapse.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_

#include "IPreSynapse.h"
#include "../synapse_handler/function/types.h" 	//SynapseType
#include <string>
#include <common/types.h> //number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" //VectorProxyBase


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PreAlphaSynapse : public IPreSynapse
{
private:
	number m_onset;

public:
	PreAlphaSynapse(const number& location, const number& onset);
	PreAlphaSynapse(const unsigned long id, const number& location, const number& onset);
	virtual ~PreAlphaSynapse();

	std::string name() const {return "PRE_ALPHA_SYNAPSE";}
	SynapseType type() const {return PRE_ALPHA_SYNAPSE;}
	number onset() const {return m_onset;}
	void set_onset(const number& onset) {m_onset = onset;}

	void update(const number& t, VectorProxyBase* up=NULL);
	bool active(const number& t, VectorProxyBase* up=NULL);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_ */
