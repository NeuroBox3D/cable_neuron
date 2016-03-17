/*
 * PreExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLITTEDSYNAPSE_HANDLER_PREEXP2SYNAPSE_H_
#define SPLITTEDSYNAPSE_HANDLER_PREEXP2SYNAPSE_H_

#include "IPreSynapse.h"											//IPreSynapse
#include "../synapse_handler/function/types.h"						//SynapseType
#include <string>													//std::string
#include <common/types.h>											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PreExp2Synapse : public IPreSynapse
{
private:
	number m_onset;

public:
	PreExp2Synapse(const number& location, const number& onset);
	PreExp2Synapse(const unsigned long id, const number& location, const number& onset);
	virtual ~PreExp2Synapse();

	std::string name() const {return "PRE_EXP2_SYNAPSE";}
	SynapseType type() const {return PRE_EXP2_SYNAPSE;}
	number onset() const {return m_onset;}
	void set_onset(const number& onset) {m_onset = onset;}

	void update(const number& t, VectorProxyBase* up=NULL);
	bool active(const number& t, VectorProxyBase* up=NULL);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_PREEXP2SYNAPSE_H_ */
