/*
 * PostAlphaSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_
#define SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_

#include <common/types.h> //number
#include "IPostSynapse.h"

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
	//ctor and dtor
	PostAlphaSynapse(const number& gMax, const number& onset, const number& tau, const number& vm, const number& e);
	virtual ~PostAlphaSynapse();

	SynapseType type() {return POST_ALPHA_SYNAPSE;}
	std::string name() {return "POST_ALPHA_SYNAPSE";}


	//setter and getter
	number gMax() const {return m_gMax;}
	void set_gMax(const number& gMax) {m_gMax = gMax;}

	number onset() const {return m_onset;}
	void set_onset(const number& onset) {m_onset = onset;}

	number tau() const {return m_tau;}
	void set_tau(const number& tau) {m_tau = tau;}

	number vm() const {return m_vm;}
	void set_vm(const number& vm) {m_vm = vm;}

	number e() const {return m_e;}
	void set_e(const number& e) {m_e = e;}


	//functionality
	number current(const number& t);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLITTEDSYNAPSE_HANDLER_POSTALPHASYNAPSE_H_ */
