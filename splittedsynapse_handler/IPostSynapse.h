/*
 * IPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_
#define SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPostSynapse {
private:

	long m_pre_synapse_id;

public:
	IPostSynapse();
	virtual ~IPostSynapse();

	void set_id(const long id) {m_pre_synapse_id = id;}
	long id() const {return m_pre_synapse_id;}

	virtual number current(const number& t) = 0;

};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_ */
