/*
 * AlphaPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ALPHA_POST_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ALPHA_POST_SYNAPSE_H

#include "post_synapse.h"   // IPostSynapse


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPostSynapse : public IPostSynapse
{
public:
    static const std::string name_string;

private:
	number m_onset;
	number m_gMax;
	number m_tau;
	number m_rev;

public:
	// ctor & dtor
	AlphaPostSynapse();		//needed for template generator
	AlphaPostSynapse(
		const number location,
		const number gMax,
		const number tau,
		const number rev);

	AlphaPostSynapse(
		const synapse_id id,
		const number location,
		const number gMax,
		const number tau,
		const number rev);

	virtual ~AlphaPostSynapse();

	// setter & getter
	void set_gMax(number gMax) {m_gMax = gMax;}
	void set_onset(number onset) {m_onset = onset;}
	void set_tau(number tau) {m_tau = tau;}
	void set_rev(number rev) {m_rev = rev;}

	number gMax() const {return m_gMax;}
	number onset() const {return m_onset;}
	number tau() const {return m_tau;}
	number rev() const {return m_rev;}

	SynapseType type() const {return ALPHA_POST_SYNAPSE;}
	const std::string& name() const {return name_string;}

	// functionality
	number current(const number& t, const number& vm);

	// serialization routines inherited from IBaseSynapse interface
    void put_to(std::ostream& os) const;			// 'put_to' == operator<<
    void get_from(std::istream& is);				// 'get_from' == operator>>

	void activate(number time);
	void deactivate();

	bool is_active(number time);

};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ALPHA_POST_SYNAPSE_H
