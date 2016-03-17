/*!
 * \file types.h
 * \brief synapse types
 *
 *  Created on: May 5, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__TYPES_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__TYPES_H__

/* \defgroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/*!
 * \brief synapses types
 */
enum SynapseType {
	EMPTY_SYNAPSE=0,
	ALPHA_SYNAPSE,
	EXP2_SYNAPSE,
	CUSTOM_SYNAPSE, //!< could be used or another for custom synapses
	ALPHA_PRE_SYNAPSE,
	ALPHA_POST_SYNAPSE,
	EXP2_PRE_SYNAPSE,
	EXP2_POST_SYNAPSE
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__FUNCTION__TYPES_H__
