/*!
 * \file types.h
 * \brief synapse types
 *
 *  Created on: May 5, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__FUNCTION__TYPES_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__FUNCTION__TYPES_H__

/* \defgroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {

/*!
 * \brief synapses types
 * @todo somehow move the association of types and integers
 *       to the synapse_dealer and automize it
 */
enum SynapseType
{
	// add new synapse types here
	// beware to keep sub-types in order:
	// all pre-synapses first, then all post-synapses
	UNDEF = 0,
	ONSET_PRE_SYNAPSE,
	THRESHOLD_PRE_SYNAPSE,
	ALPHA_POST_SYNAPSE,
	EXP2_POST_SYNAPSE
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__FUNCTION__TYPES_H__
