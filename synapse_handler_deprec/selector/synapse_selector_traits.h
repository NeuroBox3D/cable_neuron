/*!
 * \file synapse_selector_traits.h
 * \brief synapse selector traits
 *
 *  Created on: Apr 17, 2015
 *      Author: stephan
 */

/// guard
#ifndef  __H__UG__SYNAPSE_HANDLING__SYNAPSE_SELECTOR_TRAITS__
#define  __H__UG__SYNAPSE_HANDLING__SYNAPSE_SELECTOR_TRAITS__

/// includes
#include "../function/synapses.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/*!
 * \brief default impl for TSynapse
 */
template <typename TSynapse>
struct synapse_selector_trait {
	typedef TSynapse value_type;
	static number current(const value_type& v, const number& t) { return 0.; }
};


/// SYNAPSES

/*!
 * \brief NEURON's alpha synapse
 */
template <>
struct synapse_selector_trait<AlphaSynapse> {
	typedef AlphaSynapse value_type;
	static number current(const value_type& v, const number& t) {
		return v.current(t);

	}
};

/*!
 * \brief NEURON's bi-exponential synapse
 */
template <>
struct synapse_selector_trait<Exp2Syn> {
	typedef Exp2Syn value_type;
	static number current(const value_type& v, const number& t) {
		return v.current(t);
	}
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // __H__UG__SYNAPSE_HANDLING__SYNAPSE_SELECTOR_TRAITS__
