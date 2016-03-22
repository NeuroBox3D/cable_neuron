/*
 * PreExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_EXP2PRESYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_EXP2PRESYNAPSE_H_

#include "../synapse_handler/function/types.h"						//SynapseType
#include <string>													//std::string
#include <common/types.h>											//number
#include "../split_synapse_handler/IPreSynapse.h"											//IPreSynapse
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PreSynapse : public IPreSynapse
{
private:
	number m_onset;

public:
	//ctor & dtor
	Exp2PreSynapse(
			const number& location,
			const number& onset);

	Exp2PreSynapse(
			const unsigned long long id,
			const unsigned long long postsynapse_id,
			const number& location,
			const number& onset);

	virtual ~Exp2PreSynapse();

	//setter & getter
	void set_onset(const number& onset) {m_onset = onset;}
	number onset() const {return m_onset;}

	SynapseType type() const {return EXP2_PRE_SYNAPSE;}
	std::string name() const {return "EXP2_PRE_SYNAPSE";}

	//interface methods
	void update(const number& t, VectorProxyBase* up=NULL);
	bool is_active(const number& t, VectorProxyBase* up=NULL);

	//serialization interface methods
	void put_to(std::ostream& os) const;					//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2PRESYNAPSE_H_ */
