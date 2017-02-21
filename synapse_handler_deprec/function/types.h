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
 */
enum SynapseType {
	ONSET_PRE_SYNAPSE = 0,
	ALPHA_POST_SYNAPSE,
	THRESHOLD_PRE_SYNAPSE,
	EXP2_POST_SYNAPSE,
    //EMPTY_SYNAPSE,
    //ALPHA_SYNAPSE,
    //EXP2_SYNAPSE,
    //JANA_SYNAPSE_FROM_MARKUS_WITH_LOVE, // credit to M.S. for creative naming :)
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__FUNCTION__TYPES_H__
