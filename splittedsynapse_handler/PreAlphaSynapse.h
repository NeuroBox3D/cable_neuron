/*
 * PreAlphaSynapse.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_

#include "IPreSynapse.h"
#include <string>
#include <common/types.h>
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" //VectorProxyBase


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PreAlphaSynapse : public IPreSynapse
{
public:
	PreAlphaSynapse(number location);
	virtual ~PreAlphaSynapse();
	std::string name() {return "PRE_ALPHA_SYNAPSE";}
	SynapseType type() {return PRE_ALPHA_SYNAPSE;}
	void update(number t, VectorProxyBase* up);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_ */
