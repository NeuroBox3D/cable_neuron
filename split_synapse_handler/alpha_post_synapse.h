/*
 * PostAlphaSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_ALPHA_POST_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_ALPHA_POST_SYNAPSE_H_

#include <common/types.h> 							//number
#include <string>									//std::string
#include <cmath>

#include "post_synapse.h" 	//IPostSynapse
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPostSynapse : public IPostSynapse {
private:
	number m_gMax;
	number m_onset;
	number m_tau;
//	number m_vm;
	number m_e;

public:
	//ctor & dtor
	AlphaPostSynapse();		//needed for template generator
	AlphaPostSynapse(
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
//		const number& vm,
		const number& e);

	AlphaPostSynapse(
		const unsigned long long id,
		const unsigned long long presynapse_id,
		const number& location,
		const number& gMax,
		const number& onset,
		const number& tau,
//		const number& vm,
		const number& e);

	virtual ~AlphaPostSynapse();

	//setter & getter
	void set_gMax(const number& gMax) {m_gMax = gMax;}
	void set_onset(const number& onset) {m_onset = onset;}
	void set_tau(const number& tau) {m_tau = tau;}
//	void set_vm(const number& vm) {m_vm = vm;}
	void set_e(const number& e) {m_e = e;}

	number gMax() const {return m_gMax;}
	number onset() const {return m_onset;}
	number tau() const {return m_tau;}
//	number vm() const {return m_vm;}
	number e() const {return m_e;}

	SynapseType type() const {return ALPHA_POST_SYNAPSE;}
	std::string name() const {return "ALPHA_POST_SYNAPSE";}

	//post synapses are false
	bool split_type() const {return false;}

	//functionality
	number current(const number& t, const number& vm);

	//serialization from IBaseSynapse interface
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>

	void activate(number time) {m_onset = time;}
	void deactivate(number time) {m_onset = nan("");}

	bool is_active(number time);

};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_ALPHA_POST_SYNAPSE_H_ */
