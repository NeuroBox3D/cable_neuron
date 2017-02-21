/*!
 * \file synapse_selector_info_traits.h
 *
 *  Created on: Apr 17, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_SELECTOR_INFO_TRAITS_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_SELECTOR_INFO_TRAITS_H__

#include "../function/types.h"

/*! \addtogroup plugin_cable_neuron Plugin cable_neuon
 * \ingroup plugins_experimental
 * \{
 */

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

#define DECLARE_SYNAPSE_SELECTOR_INFO_TRAITS(synapseType, synapseName)\
	template <> struct synapse_selector_info_traits<synapseType> {\
		static SynapseType type_name ()	{\
			return synapseName;\
		}\
	};

template <typename TSynapse>
struct synapse_selector_info_traits {
	static SynapseType type_name ();
};

DECLARE_SYNAPSE_SELECTOR_INFO_TRAITS(Exp2Syn, EXP2_SYNAPSE);
DECLARE_SYNAPSE_SELECTOR_INFO_TRAITS(AlphaSynapse, ALPHA_SYNAPSE);

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
/// \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SELECTOR__SYNAPSE_SELECTOR_INFO_TRAITS_H__
