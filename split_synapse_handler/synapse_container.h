/*
 * synapse_container.h
 *
 *  Created on: Apr 17, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <common/types.h> 											//number
#include <cmath>

#include "alpha_post_synapse.h"
#include "alpha_pre_synapse.h"
#include "exp2_post_synapse.h"
#include "exp2_pre_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {



class ISynapseContainer {
public:
	virtual ~ISynapseContainer() {}
	virtual std::vector<IBaseSynapse*> get_synapses() = 0;
	virtual size_t size() = 0;
};



/**
 * class for parametrization of a set of AlphaSynapses, for usage in lua-scripts.
 * creates <num_synapses> pairs of a post- and presynapse, in that order. Therefore the first ID will be <start_id>
 * and the last id will be <start_id> + 2*<num_synapses> - 1
 */
class AlphaSynapses : ISynapseContainer
{
private:
	size_t m_start_id;
	number m_mean_onset;
	number m_dev_onset;
	number m_mean_gMax;
	number m_dev_gMax;
	number m_mean_tau;
	number m_dev_tau;
	number m_mean_rev;
	number m_dev_rev;

	std::vector<IBaseSynapse*> m_vSynapses;

public:
	AlphaSynapses(const size_t start_id, const size_t num_synapses)
:m_start_id(start_id),
 m_mean_onset(0),
 m_dev_onset(0),
 m_mean_gMax(0),
 m_dev_gMax(0),
 m_mean_tau(0),
 m_dev_tau(0),
 m_mean_rev(0),
 m_dev_rev(0),
 m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

	AlphaSynapses(
			const size_t start_id,
			const size_t num_synapses,
			number mean_onset,
			number dev_onset,
			number mean_gMax,
			number dev_gMax,
			number mean_tau,
			number dev_tau,
			number mean_rev,
			number dev_rev)
:m_start_id(start_id),
 m_mean_onset(mean_onset),
 m_dev_onset(dev_onset),
 m_mean_gMax(mean_gMax),
 m_dev_gMax(dev_gMax),
 m_mean_tau(mean_tau),
 m_dev_tau(dev_tau),
 m_mean_rev(mean_rev),
 m_dev_rev(dev_rev),
 m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

	void set_mean_onset(const number val) {m_mean_onset=val;}
	void set_dev_onset(const number val) {m_dev_onset=val;}
	void set_mean_gMax(const number val) {m_mean_gMax=val;}
	void set_dev_gMax(const number val) {m_dev_gMax=val;}
	void set_mean_tau(const number val) {m_mean_tau=val;}
	void set_dev_tau(const number val) {m_dev_tau=val;}
	void set_mean_rev(const number val) {m_mean_rev=val;}
	void set_dev_rev(const number val) {m_dev_rev=val;}
	size_t size() {return m_vSynapses.size();}
	size_t start_id() {return m_start_id;}

	// fixme: delete synapses in m_vSynapses!
	virtual ~AlphaSynapses() {}

	/**
	 * Assemble container
	 */
	std::vector<IBaseSynapse*> get_synapses() {

		//return m_vSynapses;

		boost::mt19937 onset_rng;
		boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

		boost::mt19937 gmax_rng;
		boost::normal_distribution<number> gmax_dist(m_mean_gMax, m_dev_gMax);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > gmax_var(gmax_rng, gmax_dist);

		boost::mt19937 tau_rng;
		boost::normal_distribution<number> tau_dist(m_mean_tau, m_dev_tau);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau_var(tau_rng, tau_dist);

		boost::mt19937 rev_rng;
		boost::normal_distribution<number> rev_dist(m_mean_rev, m_dev_rev);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > rev_var(rev_rng, rev_dist);


		//create a couple of synapses, pre and post synapse together.
		for(size_t i = 0; i < m_vSynapses.size(); i=i+2) {
			number onset = onset_var();
			number tau = tau_var();
			size_t syn_id = m_start_id + i/2;
			IBaseSynapse *s1 = new AlphaPostSynapse(syn_id, 0.0, nan(""), gmax_var(), tau, rev_var());
			IBaseSynapse *s2 = new AlphaPreSynapse(syn_id, 0.0, onset, 6 * tau);
			m_vSynapses[i] = s1;
			m_vSynapses[i+1] = s2;
		}

		return m_vSynapses;
	}
};


/**
 * class for parametrization of a set of Exp2Synapses, for usage in lua-scripts
 */
class Exp2Synapses : ISynapseContainer
{
private:
	size_t m_start_id;
	number m_mean_onset;
	number m_dev_onset;
	number m_mean_gMax;
	number m_dev_gMax;
	number m_mean_tau1;
	number m_dev_tau1;
	number m_mean_tau2;
	number m_dev_tau2;
	number m_mean_rev;
	number m_dev_rev;
	number m_mean_threshold;
	number m_dev_threshold;

	std::vector<IBaseSynapse*> m_vSynapses;

public:
	Exp2Synapses(const size_t start_id, const size_t num_synapses)
:m_start_id(start_id),
 m_mean_onset(0),
 m_dev_onset(0),
 m_mean_gMax(0),
 m_dev_gMax(0),
 m_mean_tau1(0),
 m_dev_tau1(0),
 m_mean_tau2(0),
 m_dev_tau2(0),
 m_mean_rev(0),
 m_dev_rev(0),
 m_mean_threshold(0),
 m_dev_threshold(0),
 m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

	Exp2Synapses(
			const size_t start_id,
			const size_t num_synapses,
			number mean_onset,
			number dev_onset,
			number mean_gMax,
			number dev_gMax,
			number mean_tau1,
			number dev_tau1,
			number mean_tau2,
			number dev_tau2,
			number mean_rev,
			number dev_rev,
			number mean_threshold,
			number dev_threshold
			)

:m_start_id(start_id),
 m_mean_onset(mean_onset),
 m_dev_onset(dev_onset),
 m_mean_gMax(mean_gMax),
 m_dev_gMax(dev_gMax),
 m_mean_tau1(mean_tau1),
 m_dev_tau1(dev_tau1),
 m_mean_tau2(mean_tau2),
 m_dev_tau2(dev_tau2),
 m_mean_rev(mean_rev),
 m_dev_rev(dev_rev),
 m_mean_threshold(mean_threshold),
 m_dev_threshold(dev_threshold),
 m_vSynapses(std::vector<IBaseSynapse*>(2 * num_synapses))
{}

	void set_mean_onset(const number val) {m_mean_onset=val;}
	void set_dev_onset(const number val) {m_dev_onset=val;}
	void set_mean_gMax(const number val) {m_mean_gMax=val;}
	void set_dev_gMax(const number val) {m_dev_gMax=val;}
	void set_mean_tau1(const number val) {m_mean_tau1=val;}
	void set_dev_tau1(const number val) {m_dev_tau1=val;}
	void set_mean_tau2(const number val) {m_mean_tau2=val;}
	void set_dev_tau2(const number val) {m_dev_tau2=val;}
	void set_mean_rev(const number val) {m_mean_rev=val;}
	void set_dev_rev(const number val) {m_dev_rev=val;}
	void set_mean_threshold(const number val) {m_mean_threshold=val;}
	void set_dev_threshold(const number val) {m_dev_threshold=val;}

	size_t size() {return m_vSynapses.size();}
	size_t start_id() {return m_start_id;}


    // fixme: delete synapses in m_vSynapses!
	virtual ~Exp2Synapses() {}
	/**
	 * Assembles container
	 */
	std::vector<IBaseSynapse*> get_synapses() {

		boost::mt19937 onset_rng;
		boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

		boost::mt19937 gMax_rng;
		boost::normal_distribution<number> gMax_dist(m_mean_gMax, m_dev_gMax);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > gMax_var(gMax_rng, gMax_dist);

		boost::mt19937 tau1_rng;
		boost::normal_distribution<number> tau1_dist(m_mean_tau1, m_dev_tau1);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau1_var(tau1_rng, tau1_dist);

		boost::mt19937 tau2_rng;
		boost::normal_distribution<number> tau2_dist(m_mean_tau2, m_dev_tau2);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > tau2_var(tau2_rng, tau2_dist);

		boost::mt19937 rev_rng;
		boost::normal_distribution<number> rev_dist(m_mean_rev, m_dev_rev);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > rev_var(rev_rng, rev_dist);

		boost::mt19937 threshold_rng;
		boost::normal_distribution<number> threshold_dist(m_mean_threshold, m_dev_threshold);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > threshold_var(threshold_rng, threshold_dist);

		for(size_t i = 0; i < m_vSynapses.size(); i=i+2) {
			number onset = nan("");
			number tau1 = tau1_var();
			number tau2 = tau2_var();

			number tp = (tau1*tau2)/(tau2-tau1) * std::log(tau2/tau1);
			tp = tp + 3 * tau2;
			size_t syn_id = m_start_id + i/2;
			IBaseSynapse *s1 = new Exp2PostSynapse(syn_id, 0.0, onset, gMax_var(), tau1, tau2, rev_var());
			IBaseSynapse *s2 = new Exp2PreSynapse(syn_id, 0.0, onset, tp, threshold_var());

			m_vSynapses[i] = s1;
			m_vSynapses[i+1] = s2;
		}

		return m_vSynapses;
	}
};



class AlphaSynapsePair
{
    public:
        /// constructor
        AlphaSynapsePair()
        : m_pre(new AlphaPreSynapse(0, 0.0, 0.0, 2.4e-3)),
          m_post(new AlphaPostSynapse(0, 0.0, 0.0, 1.2e-9, 4e-4, 0.0)) {};

        /// destructor
        ~AlphaSynapsePair() {delete m_pre; delete m_post;};


        /// set id
        void set_id(size_t id) {m_pre->set_id(id); m_post->set_id(id);}

        /// set onset
        void set_onset(number onset) {m_pre->set_onset(onset);}

        /// set tau (and duration, which is 6*tau)
        void set_tau(number tau) {m_post->set_tau(tau); m_pre->set_duration(6*tau);}

        /// set maximal conductance
        void set_gMax(number gMax) {m_post->set_gMax(gMax);}

        /// set reversal potential
        void set_reversal_potential(number ve) {m_post->set_rev(ve);}

        /// get pre-synapse
        AlphaPreSynapse* pre_synapse() {return m_pre;}

        /// get post-synapse
        AlphaPostSynapse* post_synapse() {return m_post;}

    private:
        AlphaPreSynapse* m_pre;
        AlphaPostSynapse* m_post;
};



} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSE_CONTAINER_H_ */
