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
#include "base_synapse.h"											//IBaseSynapse
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPreSynapse : public IBaseSynapse
{
	SYNAPSE_ID m_id;					//own presynapse id / alternative: shared id if there is a 1:1 relation?
	SYNAPSE_ID m_postsynapse_id;		//postsynapse id

public:
	//ctor & dtor
	IPreSynapse();						//needed for SynapseDealer
	IPreSynapse(
			const SYNAPSE_ID id,
			const SYNAPSE_ID postsynapse_id,
			const number& location);

	virtual ~IPreSynapse();

	//setter & getter
	void set_id(const SYNAPSE_ID id) {m_id = id;}
	void set_postsynapse_id(const SYNAPSE_ID id) {m_postsynapse_id = id;}

	SYNAPSE_ID id() const {return m_id;}
	SYNAPSE_ID postsynapse_id() const {return m_postsynapse_id;}

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

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_ */
