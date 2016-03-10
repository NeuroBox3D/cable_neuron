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
	long m_id;
public:
	IPreSynapse(const number& location);
	virtual ~IPreSynapse();

	number location() const {return m_location;}
	void set_location(const number& location) {m_location = location;}

	long id() const {return m_id;}
	void set_id(const int id) {m_id = id;}

	virtual std::string name() const = 0;
	virtual SynapseType type() const = 0;
	virtual void update(const number& t, VectorProxyBase* up=NULL) = 0;
	virtual bool active(const number& t, VectorProxyBase* up=NULL) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_ */
