/*
 * IPreSynapse.h
 *
 *  Created on: Feb 29, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__PRE_SYNAPSE_H_
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__PRE_SYNAPSE_H_

#include "base_synapse.h"   // IBaseSynapse, synapse_id
#include <vector>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPreSynapse : public IBaseSynapse
{


public:
	//ctor & dtor
	IPreSynapse();						//needed for SynapseDealer
	IPreSynapse(
			const synapse_id id,
			const number location);

	virtual ~IPreSynapse();

	//setter & getter
	virtual void update(const number& t, const std::vector<number>& u) = 0;
	virtual bool is_active(const number& t) = 0;

	virtual bool is_presynapse() const {return true;}
	virtual bool is_postsynapse() const {return false;}

	//from serialization interface IBaseSynapse
	virtual void put_to(std::ostream& os) const = 0;			//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;

//	virtual bool fire(number time, SYNAPSE_ID& post_syn_id) = 0;
//	virtual bool cooldown(number time, SYNAPSE_ID& post_syn_id) = 0;
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__PRE_SYNAPSE_H_
