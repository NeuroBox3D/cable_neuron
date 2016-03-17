/*
 * PreAlphaSynapse.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_

#include "../synapse_handler/function/types.h" 						//SynapseType
#include <string>													//std::string
#include <common/types.h> 											//number
#include "../split_synapse_handler/IPreSynapse.h"					//IPreSynapse
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PreAlphaSynapse : public IPreSynapse
{
private:
	number m_onset;

public:
	//ctor & dtor
	PreAlphaSynapse(
			const number& location,
			const number& onset);

	PreAlphaSynapse(
			const unsigned long id,
			const unsigned long postsynapse_id,
			const number& location,
			const number& onset);

	virtual ~PreAlphaSynapse();

	//setter & getter
	SynapseType type() const {return ALPHA_PRE_SYNAPSE;}
	std::string name() const {return "ALPHA_PRE_SYNAPSE";}

	void set_onset(const number& onset) {m_onset = onset;}
	number onset() const {return m_onset;}

	//Interface dummy's
	void update(const number& t, VectorProxyBase* up=NULL);
	bool is_active(const number& t, VectorProxyBase* up=NULL);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_ */
