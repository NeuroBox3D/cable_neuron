/*
 * PostExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__EXP2_POST_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__EXP2_POST_SYNAPSE_H

#include "post_synapse.h"   // IPostSynapse


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PostSynapse : public IPostSynapse
{
    public:
        static const std::string name_string;

private:
	number m_onset; ///<onset (s)
	number m_gMax;  ///< conductance (S)
	number m_tau1;  ///< time constant (s)
	number m_tau2;  ///< time constant (s)
	number m_rev;   ///< reversal potential (V)


public:
	//ctor & dtor
	Exp2PostSynapse();  // needed for template generator
	Exp2PostSynapse(
			const number location,
			const number gMax,
			const number tau1,
			const number tau2,
			const number rev
			);

	Exp2PostSynapse(
			const synapse_id id,
			const number location,
			const number gMax,
			const number tau1,
			const number tau2,
			const number rev
			);

	virtual ~Exp2PostSynapse();

	//setter & getter
	void set_onset(number onset) {m_onset=onset;}
	void set_gMax(number gmax) {m_gMax=gmax;}
	void set_tau1(number tau1) {m_tau1=tau1;}
	void set_tau2(number tau2) {m_tau2=tau2;}
	void set_rev(number rev) {m_rev=rev;}


	number onset() const {return m_onset;}
	number gMax() const {return m_gMax;}
	number tau1() const {return m_tau1;}
	number tau2() const {return m_tau2;}
	number rev() const {return m_rev;}


	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	const std::string& name() const {return name_string;}

	//functionality
	number current(const number& t, const number& vm);

	//serialization interface methods
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>

	void activate(number time);
	void deactivate();

	bool is_active(number time);
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__EXP2_POST_SYNAPSE_H
