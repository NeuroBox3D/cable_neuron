/*
 * IBaseSynapse.h
 *
 *  Created on: Mar 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_IBASESYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_IBASESYNAPSE_H_

#include <iostream>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IBaseSynapse {

public:
	IBaseSynapse(){}
	virtual ~IBaseSynapse(){}

	//todo: maybe protected for polymorphic ops and private for global friend ops?
	friend std::ostream& operator<<(std::ostream& os, const IBaseSynapse* s);
	friend std::istream& operator>>(std::istream& is, IBaseSynapse* s);
	virtual void put_to(std::ostream& os) const = 0;							//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;								//'get_from' == operator>>
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_IBASESYNAPSE_H_ */
