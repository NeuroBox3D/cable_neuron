/*
 * PreExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_EXP2_PRE_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_EXP2_PRE_SYNAPSE_H_

#include "../synapse_handler/function/types.h"						//SynapseType
#include <string>													//std::string
#include <common/types.h>											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase
#include "pre_synapse.h"											//IPreSynapse
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/normal_distribution.hpp>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PreSynapse : public IPreSynapse
{
private:
	number m_onset;
	number m_duration;
	number m_threshold;

public:
	//ctor & dtor
	Exp2PreSynapse();					//needed for template generator
	Exp2PreSynapse(
			const number& location,
			const number& onset,
			const number& duration
			);

	Exp2PreSynapse(
			const SYNAPSE_ID id,
			const SYNAPSE_ID postsynapse_id,
			const number& location,
			const number& onset,
			const number& duration
			);

	virtual ~Exp2PreSynapse();

	//setter & getter
	void set_onset(const number& onset) {m_onset = onset;}
	void set_duration(const number& duration) {m_duration = duration;}

	number onset() const {return m_onset;}
	number duration() const {return m_duration;}

	SynapseType type() const {return EXP2_PRE_SYNAPSE;}
	std::string name() const {return "EXP2_PRE_SYNAPSE";}

	//pre synapses are true
	bool split_type() const {return true;}

	//interface methods
	void update(const number& t, const std::vector<number>& u);
	bool is_active(const number& t);

	//serialization interface methods
	void put_to(std::ostream& os) const;					//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>

//	bool fire(number time, SYNAPSE_ID& post_syn_id);
//	bool cooldown(number time, SYNAPSE_ID& post_syn_id);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2_PRE_SYNAPSE_H_ */
