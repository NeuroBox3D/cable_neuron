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
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase
#include "pre_synapse.h"					//IPreSynapse


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPreSynapse : public IPreSynapse
{
private:
	number m_onset;

public:
	//ctor & dtor
	AlphaPreSynapse() {} 					//needed for template generator
	AlphaPreSynapse(
			const number& location,
			const number& onset);

	AlphaPreSynapse(
			const unsigned long long id,
			const unsigned long long postsynapse_id,
			const number& location,
			const number& onset);

	virtual ~AlphaPreSynapse();

	//setter & getter
	SynapseType type() const {return ALPHA_PRE_SYNAPSE;}
	std::string name() const {return "ALPHA_PRE_SYNAPSE";}

	void set_onset(const number& onset) {m_onset = onset;}
	number onset() const {return m_onset;}

	//Interface dummy's
	void update(const number& t, VectorProxyBase* up=NULL);
	bool is_active(const number& t, VectorProxyBase* up=NULL);

	//serialization from IBaseSynapse interface
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_ */
