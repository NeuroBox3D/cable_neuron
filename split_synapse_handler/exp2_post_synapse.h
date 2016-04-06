/*
 * PostExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_

#include <string>													//std::string
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <common/types.h>											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase
#include "post_synapse.h"											//IPostSynapse

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PostSynapse : public IPostSynapse
{
private:
	number m_tau1;	///< time constant (s)
	number m_tau2;	///< time constant (s)
	number m_e;		///< reversal potential (V)
	number m_w;		///< conductance (S)
	number m_vm;	///< membrane potential (V)

public:
	//ctor & dtor
	Exp2PostSynapse();					//needed for template generator
	Exp2PostSynapse(
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w,
			const number& vm);

	Exp2PostSynapse(
			const unsigned long long id,
			const unsigned long long presynapse_id,
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w,
			const number& vm);

	virtual ~Exp2PostSynapse();

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

	void set_activation_timing(std::vector<number> timings);

	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	std::string name() const {return "EXP2_POST_SYNAPSE";}

	//post synapses are false
	bool split_type() const {return false;}

	//functionality
	number current(const number& t);

	//serialization interface methods
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_ */
