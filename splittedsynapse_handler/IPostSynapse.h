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
	number m_gMax;
	number m_tau;
	number m_vm;
	number m_e;

	int m_pre_synapse_id;

public:
	IPostSynapse(const number& gMax,const number& tau,const number& vm,const number& e);
	virtual ~IPostSynapse();

	virtual number current(const number& t);

};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_IPOSTSYNAPSE_H_ */
