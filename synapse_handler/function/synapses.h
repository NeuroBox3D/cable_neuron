/*!
 * \file synapse.h
 * \brief synapses
 *
 *  Created on: Jan 19, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSES_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSES_H__

/// includes
#include <map>
#include <cmath>

#include <common/util/smart_pointer.h>

#include "../synapse_factory.h"
#include "../synapse.h"

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

	/*!
	 * \brief alpha-synapse, i. e. input synapse / unilateral (nf) synapse
	 */
	class AlphaSynapse : public ISynapse {
	private:
		number m_gMax;
		number m_onset;
		number m_tau;
		number m_e;
		number m_vm;

	public:
		number current(const number& t) const {
			if (t < m_onset) return 0.0;
			return m_gMax * (t - m_onset)/m_tau * std::exp(-(t - m_onset - m_tau)/m_tau)
				  * (m_vm - m_e);
		}

		void set_vm(const number& vm) {
			m_vm = vm;
		}

		void set_parameters(number gmax, number onset, number tau, number e) {
			m_gMax = gmax;
			m_onset = onset;
			m_tau = tau;
			m_e = e;
		}

		AlphaSynapse() : m_gMax(0), m_onset(0), m_tau(1), m_e(0), m_vm(0) {

		}

		AlphaSynapse(number gmax, number onset, number tau, number e, number vm) :
			m_gMax(gmax), m_onset(onset), m_tau(tau), m_e(e), m_vm(vm) {

		}

		virtual ~AlphaSynapse() {

		}
	};

	/*!
	 * \brief bi-exponential synapse, i. e. functional synapse / bilateral synapse
	 */
	class Exp2Syn : public ISynapse {
	private:
		number m_tau1;
		number m_tau2;
		number m_e;
		number m_w;
		number m_vm;

	public:
		number current(const number& t) const {
			number tp = (m_tau1*m_tau2)/(m_tau2 - m_tau1) * std::log(m_tau2/m_tau1);	// time of maximal current
			number factor = 1.0 / (std::exp(-tp/m_tau2) - std::exp(-tp/m_tau1));		// normalization factor
			number i = m_w * factor * (m_vm - m_e) * (std::exp(-t/m_tau2) - std::exp(-t/m_tau1));
			return i; //!< i: current [nA]
		}

		void set_vm(const number& vm) {
			m_vm = vm;
		}

		void set_parameters(number tau1, number tau2, number e, number w) {
			m_tau1 = tau1;
			m_tau2 = tau2;
			m_e = e;
			m_w = w;
		}

		Exp2Syn() : m_tau1(1), m_tau2(1), m_e(0), m_w(0), m_vm(0) {

		}

		Exp2Syn(number tau1, number tau2, number e, number w, number vm) : m_tau1(tau1), m_tau2(tau2), m_e(e), m_w(w), m_vm(vm) {

		}

		virtual ~Exp2Syn() {

		}
	};


// ///////////////////////////////
// factories 					//
// ///////////////////////////////
	/*!
	 * \brief a synapse factory for alpha synapses
	 */
	class AlphaSynapseFactory : public ISynapseFactory {
	public:
		AlphaSynapse* create_synapse() {
			return new AlphaSynapse();
		}

		AlphaSynapseFactory() {

		}
		virtual ~AlphaSynapseFactory() {

		}
	};

	/*!
	 * \brief a synapse factory for bi-exponential synapses
	 */
	class Exp2SynapseFactory : public ISynapseFactory {
	public:
		Exp2Syn* create_synapse(number tau1, number tau2, number e, number w, number vm) {
			return new Exp2Syn(tau1, tau2, e, w, vm);
		}

		Exp2Syn* create_synapse() {
			return new Exp2Syn();
		}

		virtual ~Exp2SynapseFactory() {
		}
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__SYNAPSES_H__
