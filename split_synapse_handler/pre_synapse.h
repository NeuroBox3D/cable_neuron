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
	unsigned long long m_id;					//own presynapse id / alternative: shared id if there is a 1:1 relation?
	unsigned long long m_postsynapse_id;		//postsynapse id
	number m_location;					//location on edge

public:
	//ctor & dtor
	IPreSynapse() {}						//needed for SynapseDealer
	IPreSynapse(
			const unsigned long long id,
			const unsigned long long postsynapse_id,
			const number& location);

	virtual ~IPreSynapse();

	//setter & getter
	void set_id(const unsigned long long id) {m_id = id;}
	void set_postsynapse_id(const unsigned long long id) {m_postsynapse_id = id;}
	void set_location(const number& location) {m_location = location;}

	unsigned long long id() const {return m_id;}
	unsigned long long postsynapse_id() const {return m_postsynapse_id;}
	number location() const {return m_location;}

	//virtual interface methods
	/**
	 * important: pre synapse implementations have to end with '_PRE_SYNAPSE'
	 * e.g. EXP2_PRE_SYNAPSE
	 */
	virtual std::string name() const = 0;
	virtual SynapseType type() const = 0;

	virtual void update(const number& t, VectorProxyBase* up=NULL) = 0;
	virtual bool is_active(const number& t, VectorProxyBase* up=NULL) = 0;

	//from serialization interface IBaseSynapse
	virtual void put_to(std::ostream& os) const = 0;			//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_IPRESYNAPSE_H_ */
