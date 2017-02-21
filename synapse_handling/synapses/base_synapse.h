/*
 * IBaseSynapse.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__BASE_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__BASE_SYNAPSE_H

#include "../types.h"
#include "common/types.h"
#include <iostream>
#include <string>

typedef unsigned long long synapse_id;

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IBaseSynapse {
private:
	synapse_id m_id;							//own presynapse id / alternative: shared id if there is a 1:1 relation?
	number m_location;							//location on edge

public:

	IBaseSynapse()
	:m_id(0),
	 m_location(0)
	{}

	IBaseSynapse(synapse_id id, number location)
	:m_id(id),
	 m_location(location)
	{}

	virtual ~IBaseSynapse(){}

	//todo: maybe protected for polymorphic ops and private for global friend ops?
	friend std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s);
	friend std::istream& operator>>(std::istream& is, IBaseSynapse* s);
	virtual void put_to(std::ostream& os) const = 0;							//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;								//'get_from' == operator>>

	virtual const std::string& name() const = 0;
	virtual SynapseType type() const = 0;
	virtual bool is_presynapse() const = 0;
	virtual bool is_postsynapse() const = 0;

	void set_location(const number& loc) {m_location = loc;}
	number location() const {return m_location;}

	void set_id(const synapse_id id) {m_id = id;}
	synapse_id id() const {return m_id;}
};


} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__BASE_SYNAPSE_H
