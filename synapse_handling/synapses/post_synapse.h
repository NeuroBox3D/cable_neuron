/*
 * IPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__POST_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__POST_SYNAPSE_H

#include "common/types.h"   // number etc.
#include "base_synapse.h"   // IBaseSynapse, synapse_id


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPostSynapse : public IBaseSynapse
{
public:
	//ctor & dtor
	IPostSynapse(); 							//needed for SynapseDealer
	IPostSynapse(
			const synapse_id id,
			const number location
	);

	virtual ~IPostSynapse();

	//setter & getter

	virtual bool is_presynapse() const {return false;}
	virtual bool is_postsynapse() const {return true;}

	/// (outward) current of a synapse
	virtual number current(const number& t, const number& vm) const = 0;

	//from serialization interface IBaseSynapse
	virtual void put_to(std::ostream& os) const = 0;			//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;				//'get_from' == operator>>

	virtual void activate(number time) = 0;
	virtual void deactivate() = 0;
	virtual bool is_active(number time) = 0;
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__POST_SYNAPSE_H
