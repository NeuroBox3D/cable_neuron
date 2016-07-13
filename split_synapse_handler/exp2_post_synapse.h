/*
 * PostExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_

#include <string>													//std::string
#include <cmath>
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <common/types.h>											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	//VectorProxyBase
#include "post_synapse.h"											//IPostSynapse
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/normal_distribution.hpp>

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class Exp2PostSynapse : public IPostSynapse
{
private:
	number m_onset; ///onset (s)
	number m_gMax;		///< conductance (S)
	number m_tau1;	///< time constant (s)
	number m_tau2;	///< time constant (s)
	number m_rev;		///< reversal potential (V)



public:
	//ctor & dtor
	Exp2PostSynapse();					//needed for template generator
	Exp2PostSynapse(
			const number& location,
			const number& onset,
			const number& gMax,
			const number& tau1,
			const number& tau2,
			const number& rev
			);

	Exp2PostSynapse(
			const SYNAPSE_ID id,
			const SYNAPSE_ID presynapse_id,
			const number& location,
			const number& onset,
			const number& gMax,
			const number& tau1,
			const number& tau2,
			const number& rev
			);

	virtual ~Exp2PostSynapse();

	//setter & getter
	void set_onset(const number& onset) {m_onset=onset;}
	void set_gMax(const number& gmax) {m_gMax=gmax;}
	void set_tau1(const number& tau1) {m_tau1=tau1;}
	void set_tau2(const number& tau2) {m_tau2=tau2;}
	void set_rev(const number& rev) {m_rev=rev;}


	number onset() const {return m_onset;}
	number gMax() const {return m_gMax;}
	number tau1() const {return m_tau1;}
	number tau2() const {return m_tau2;}
	number rev() const {return m_rev;}


	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	std::string name() const {return "EXP2_POST_SYNAPSE";}

	//functionality
	number current(const number& t, const number& vm);

	//serialization interface methods
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>

	void activate(number time) {m_onset = time;}
	void deactivate(number time) {m_onset = nan("");}

	bool is_active(number time);
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_ */
