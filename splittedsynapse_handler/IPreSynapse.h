/*
 * IPreSynapse.h
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_

#include <string>								//std::string
#include "../synapse_handler/function/types.h" 	//SynapseType
#include <common/types.h>						//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" //VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPreSynapse {
	number m_location;
public:
	IPreSynapse(number location);
	number location() {return m_location;}

	virtual ~IPreSynapse();
	virtual std::string name() = 0;
	virtual SynapseType type() = 0;
	virtual void update(number t, VectorProxyBase* up) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_ */
