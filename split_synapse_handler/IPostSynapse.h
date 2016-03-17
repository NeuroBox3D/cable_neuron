/*
 * IPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_
#define SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_

#include <common/types.h> //number
#include "../synapse_handler/function/types.h" 	//SynapseType
#include <string>
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" //VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPostSynapse {
private:
	unsigned long m_presynapse_id;
	number m_location;

public:
	//ctor & dtor
	IPostSynapse(const number& location);
	IPostSynapse(const unsigned long presynapseid, const number& location);
	virtual ~IPostSynapse();

	//setter & getter
	void set_presynapse_id(const unsigned long id) {m_presynapse_id = id;}
	void set_location(const number& loc) {m_location = loc;}
	unsigned long presynapse_id() const {return m_presynapse_id;}
	number location() const {return m_location;}

	//virtual interface methods
	virtual SynapseType type() const = 0;
	virtual std::string name() const = 0;

	virtual number current(const number& t) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_ */
