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

#include "post_synapse.h" 	//IPostSynapse
#include "synapse_container.h"
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

};


/**
 * class for parametrization of a set of AlphaSynapses, for usage in lua-scripts
 */
class AlphaPostSynapses : ISynapseContainer
{
private:
	number m_mean_gMax;
	number m_dev_gMax;
	number m_mean_onset;
	number m_dev_onset;
	number m_mean_tau;
	number m_dev_tau;
	number m_mean_e;
	number m_dev_e;

	std::vector<IBaseSynapse*> m_vSynapses;

	bool m_parametrizised;

public:
	AlphaPostSynapses(const size_t num_synapses)
:m_mean_gMax(0),
 m_dev_gMax(0),
 m_mean_onset(0),
 m_dev_onset(0),
 m_mean_tau(0),
 m_dev_tau(0),
 m_mean_e(0),
 m_dev_e(0),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	AlphaPostSynapses(
			const size_t num_synapses,
			number mean_gMax,
			number dev_gMax,
			number mean_onset,
			number dev_onset,
			number mean_tau,
			number dev_tau,
			number mean_e,
			number dev_e)
:m_mean_gMax(mean_gMax),
 m_dev_gMax(dev_gMax),
 m_mean_onset(mean_onset),
 m_dev_onset(dev_onset),
 m_mean_tau(mean_tau),
 m_dev_tau(dev_tau),
 m_mean_e(mean_e),
 m_dev_e(dev_e),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	void set_mean_gMax(const number val) {m_mean_gMax=val;}
	void set_dev_gMax(const number val) {m_dev_gMax=val;}
	void set_mean_onset(const number val) {m_mean_onset=val;}
	void set_dev_onset(const number val) {m_dev_onset=val;}
	void set_mean_tau(const number val) {m_mean_tau=val;}
	void set_dev_tau(const number val) {m_dev_tau=val;}
	void set_mean_e(const number val) {m_mean_e=val;}
	void set_dev_e(const number val) {m_dev_e=val;}
	size_t size() {return m_vSynapses.size();}

	virtual ~AlphaPostSynapses() {}
	/**
	 * Will be called by C++ functions
	 */
	std::vector<IBaseSynapse*> get_synapses() {



		//return m_vSynapses;

		if(!m_parametrizised)
			UG_THROW("No parameters for AlphaPostSynapses");

		boost::mt19937 gmax_rng;
		boost::normal_distribution<number> gmax_dist(m_mean_gMax, m_dev_gMax);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > gmax_var(gmax_rng, gmax_dist);

		boost::mt19937 onset_rng;
		boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

		boost::mt19937 tau_rng;
		boost::normal_distribution<number> tau_dist(m_mean_tau, m_dev_tau);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau_var(tau_rng, tau_dist);

		boost::mt19937 e_rng;
		boost::normal_distribution<number> e_dist(m_mean_e, m_dev_e);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > e_var(e_rng, e_dist);


		for(size_t i = 0; i < m_vSynapses.size(); ++i) {
			IBaseSynapse *s = new AlphaPostSynapse(0.0, gmax_var(), onset_var(), tau_var(), e_var());
			m_vSynapses[i] = s;
		}

		return m_vSynapses;
	}
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_ALPHA_POST_SYNAPSE_H_ */
