/*
 * PostExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_EXP2POSTSYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_EXP2POSTSYNAPSE_H_

#include <string>													//std::string
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <common/types.h>											//number
#include "../split_synapse_handler/IPostSynapse.h"											//IPostSynapse
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class PostExp2Synapse : public IPostSynapse
{
private:
	number m_tau1;	///< time constant (s)
	number m_tau2;	///< time constant (s)
	number m_e;		///< reversal potential (V)
	number m_w;		///< conductance (S)
	number m_vm;	///< membrane potential (V)

public:
	//ctor & dtor
	PostExp2Synapse(
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w,
			const number& vm);

	PostExp2Synapse(
			const unsigned long id,
			const unsigned long presynapse_id,
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w,
			const number& vm);

	virtual ~PostExp2Synapse();

	//setter & getter
	void set_tau1(const number& tau1) {m_tau1=tau1;}
	void set_tau2(const number& tau2) {m_tau2=tau2;}
	void set_e(const number& e) {m_e=e;}
	void set_w(const number& w) {m_w=w;}
	void set_vm(const number& vm) {m_vm=vm;}

	number tau1() const {return m_tau1;}
	number tau2() const {return m_tau2;}
	number e() const {return m_e;}
	number w() const {return m_w;}
	number vm() const {return m_vm;}

	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	std::string name() const {return "EXP2_POST_SYNAPSE";}

	//functionality
	number current(const number& t);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2POSTSYNAPSE_H_ */
