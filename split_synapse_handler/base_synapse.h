/*
 * IBaseSynapse.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_BASE_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_BASE_SYNAPSE_H_

#include <iostream>
#include "../synapse_handler/function/types.h"
#include "common/types.h"
#include <vector>

typedef unsigned long long SYNAPSE_ID;

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IBaseSynapse {
private:
	SYNAPSE_ID m_id;							//own presynapse id / alternative: shared id if there is a 1:1 relation?
	number m_location;							//location on edge

public:

	IBaseSynapse()
	:m_id(0),
	 m_location(0)
	{}

	IBaseSynapse(SYNAPSE_ID id, number location)
	:m_id(id),
	 m_location(location)
	{}

	virtual ~IBaseSynapse(){}

	//todo: maybe protected for polymorphic ops and private for global friend ops?
	friend std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s);
	friend std::istream& operator>>(std::istream& is, IBaseSynapse* s);
	virtual void put_to(std::ostream& os) const = 0;							//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;								//'get_from' == operator>>

	virtual std::string name() const = 0;
	virtual SynapseType type() const = 0;
	virtual bool is_presynapse() const = 0;
	virtual bool is_postsynapse() const = 0;

	void set_location(const number& loc) {m_location = loc;}
	number location() const {return m_location;}

	void set_id(const SYNAPSE_ID id) {m_id = id;}
	SYNAPSE_ID id() const {return m_id;}
};


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_BASE_SYNAPSE_H_ */
