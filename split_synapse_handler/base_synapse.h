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

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IBaseSynapse {
private:
	number m_location;							//location on edge

public:
	IBaseSynapse()
	:m_location(0)
	{}

	IBaseSynapse(number location)
	:m_location(location)
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

	void set_location(const number& loc) {m_location = loc;}
	number location() const {return m_location;}
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_BASE_SYNAPSE_H_ */
