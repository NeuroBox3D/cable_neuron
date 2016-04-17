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
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "synapse_container.h"

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
//	number m_vm;	///< membrane potential (V)

public:
	//ctor & dtor
	Exp2PostSynapse();					//needed for template generator
	Exp2PostSynapse(
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w
//			const number& vm
			);

	Exp2PostSynapse(
			const unsigned long long id,
			const unsigned long long presynapse_id,
			const number& location,
			const number& tau1,
			const number& tau2,
			const number& e,
			const number& w
//			const number& vm
			);

	virtual ~Exp2PostSynapse();

	//setter & getter
	void set_tau1(const number& tau1) {m_tau1=tau1;}
	void set_tau2(const number& tau2) {m_tau2=tau2;}
	void set_e(const number& e) {m_e=e;}
	void set_w(const number& w) {m_w=w;}
//	void set_vm(const number& vm) {m_vm=vm;}

	number tau1() const {return m_tau1;}
	number tau2() const {return m_tau2;}
	number e() const {return m_e;}
	number w() const {return m_w;}
//	number vm() const {return m_vm;}

	SynapseType type() const {return EXP2_POST_SYNAPSE;}
	std::string name() const {return "EXP2_POST_SYNAPSE";}

	//post synapses are false
	bool split_type() const {return false;}

	//functionality
	number current(const number& t, const number& vm);

	//serialization interface methods
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};


/**
 * class for parametrization of a set of Exp2PostSynapses, for usage in lua-scripts
 */
class Exp2PostSynapses : ISynapseContainer
{
private:
	number m_mean_tau1;
	number m_dev_tau1;
	number m_mean_tau2;
	number m_dev_tau2;
	number m_mean_e;
	number m_dev_e;
	number m_mean_w;
	number m_dev_w;

	std::vector<IBaseSynapse*> m_vSynapses;

	bool m_parametrizised;

public:
	Exp2PostSynapses(const size_t num_synapses)
:m_mean_tau1(0),
 m_dev_tau1(0),
 m_mean_tau2(0),
 m_dev_tau2(0),
 m_mean_e(0),
 m_dev_e(0),
 m_mean_w(0),
 m_dev_w(0),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	Exp2PostSynapses(
			const size_t num_synapses,
			number mean_tau1,
			number dev_tau1,
			number mean_tau2,
			number dev_tau2,
			number mean_e,
			number dev_e,
			number mean_w,
			number dev_w)

:m_mean_tau1(mean_tau1),
 m_dev_tau1(dev_tau1),
 m_mean_tau2(mean_tau2),
 m_dev_tau2(dev_tau2),
 m_mean_e(mean_e),
 m_dev_e(dev_e),
 m_mean_w(mean_w),
 m_dev_w(dev_w),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	void set_mean_tau1(const number val) {m_mean_tau1=val;}
	void set_dev_tau1(const number val) {m_dev_tau1=val;}
	void set_mean_tau2(const number val) {m_mean_tau2=val;}
	void set_dev_tau2(const number val) {m_dev_tau2=val;}
	void set_mean_e(const number val) {m_mean_e=val;}
	void set_dev_e(const number val) {m_dev_e=val;}
	void set_mean_w(const number val) {m_mean_w=val;}
	void set_dev_w(const number val) {m_dev_w=val;}

	size_t size() {return m_vSynapses.size();}

	virtual ~Exp2PostSynapses() {}
	/**
	 * Will be called by C++ functions
	 */
	std::vector<IBaseSynapse*> get_synapses() {



		//return m_vSynapses;

		if(!m_parametrizised)
			UG_THROW("No parameters for AlphaPostSynapses");

		boost::mt19937 tau1_rng;
		boost::normal_distribution<number> tau1_dist(m_mean_tau1, m_dev_tau1);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau1_var(tau1_rng, tau1_dist);

		boost::mt19937 tau2_rng;
		boost::normal_distribution<number> tau2_dist(m_mean_tau2, m_dev_tau2);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau2_var(tau2_rng, tau2_dist);

		boost::mt19937 w_rng;
		boost::normal_distribution<number> w_dist(m_mean_w, m_dev_w);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > w_var(w_rng, w_dist);

		boost::mt19937 e_rng;
		boost::normal_distribution<number> e_dist(m_mean_e, m_dev_e);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > e_var(e_rng, e_dist);


		for(size_t i = 0; i < m_vSynapses.size(); ++i) {
			IBaseSynapse *s = new Exp2PostSynapse(0.0, tau1_var(), tau2_var(), e_var(), w_var());
			m_vSynapses[i] = s;
		}

		return m_vSynapses;
	}
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_EXP2_POST_SYNAPSE_H_ */
