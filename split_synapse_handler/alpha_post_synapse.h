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
#include <boost/lexical_cast.hpp>
#include <boost/random/normal_distribution.hpp>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPostSynapse : public IPostSynapse
{
private:
	number m_onset;
	number m_gMax;
	number m_tau;
	number m_rev;

public:
	//ctor & dtor
	AlphaPostSynapse();		//needed for template generator
	AlphaPostSynapse(
		const number location,
		const number onset,
		const number gMax,
		const number tau,
		const number rev);

	AlphaPostSynapse(
		const SYNAPSE_ID id,
		const number location,
		const number onset,
		const number gMax,
		const number tau,
		const number rev);

	virtual ~AlphaPostSynapse();

	//setter & getter
	void set_gMax(number gMax) {m_gMax = gMax;}
	void set_onset(number onset) {m_onset = onset;}
	void set_tau(number tau) {m_tau = tau;}
	void set_rev(number rev) {m_rev = rev;}

	number gMax() const {return m_gMax;}
	number onset() const {return m_onset;}
	number tau() const {return m_tau;}
	number rev() const {return m_rev;}

	SynapseType type() const {return ALPHA_POST_SYNAPSE;}
	std::string name() const {return "ALPHA_POST_SYNAPSE";}

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
