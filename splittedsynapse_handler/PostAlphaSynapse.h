/*
 * PostAlphaSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_
#define SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_

#include <common/types.h> //number

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PostAlphaSynapse : public IPostSynapse {
private:
	number m_gMax;
	number m_onset;
	number m_tau;
	number m_vm;
	number m_e;

public:
	PostAlphaSynapse(const number& gMax, const number& onset, const number& tau, const number& tau, const number& e);
	virtual ~PostAlphaSynapse();

	number current(const number& t);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_ */
