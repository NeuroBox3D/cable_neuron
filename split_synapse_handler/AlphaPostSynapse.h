/*
 * PostAlphaSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_ALPHAPOSTSYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_ALPHAPOSTSYNAPSE_H_

#include <common/types.h> 							//number
#include <string>									//std::string
#include "../split_synapse_handler/IPostSynapse.h" 	//IPostSynapse

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPostSynapse : public IPostSynapse {
private:
	number m_gMax;
	number m_onset;
	number m_tau;
	number m_vm;
	number m_e;

public:
	//ctor & dtor
	AlphaPostSynapse(
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
		const number& vm,
		const number& e);

	AlphaPostSynapse(
		const unsigned long id,
		const unsigned long presynapse_id,
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
		const number& vm,
		const number& e);

	virtual ~AlphaPostSynapse();

	//setter & getter
	void set_gMax(const number& gMax) {m_gMax = gMax;}
	void set_onset(const number& onset) {m_onset = onset;}
	void set_tau(const number& tau) {m_tau = tau;}
	void set_vm(const number& vm) {m_vm = vm;}
	void set_e(const number& e) {m_e = e;}

	number gMax() const {return m_gMax;}
	number onset() const {return m_onset;}
	number tau() const {return m_tau;}
	number vm() const {return m_vm;}
	number e() const {return m_e;}

	SynapseType type() const {return ALPHA_POST_SYNAPSE;}
	std::string name() const {return "ALPHA_POST_SYNAPSE";}

	//functionality
	number current(const number& t);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_ALPHAPOSTSYNAPSE_H_ */
