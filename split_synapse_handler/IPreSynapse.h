/*
 * IPreSynapse.h
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_

#include <string>													//std::string
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <common/types.h>											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPreSynapse
{
	unsigned long m_id;
	unsigned long m_postsynapse_id;
	number m_location;

public:
	//ctor & dtor
	IPreSynapse();
	IPreSynapse(const number& location);
	IPreSynapse(
			const unsigned long id,
			const unsigned long postsynapse_id,
			const number& location);

	virtual ~IPreSynapse();

	//setter & getter
	void set_id(const unsigned long id) {m_id = id;}
	void set_postsynapse_id(const unsigned long id) {m_postsynapse_id = id;}
	void set_location(const number& location) {m_location = location;}

	unsigned long id() const {return m_id;}
	unsigned long postsynapse_id() const {return m_postsynapse_id;}
	number location() const {return m_location;}

	//virtual interface methods
	virtual std::string name() const = 0;
	virtual SynapseType type() const = 0;

	virtual void update(const number& t, VectorProxyBase* up=NULL) = 0;
	virtual bool active(const number& t, VectorProxyBase* up=NULL) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_ */
