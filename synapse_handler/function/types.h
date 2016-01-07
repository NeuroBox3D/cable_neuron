/*!
 * \file types.h
 * \brief synapse types
 *
 *  Created on: May 5, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__SYNAPSE_HANDLER__TYPES__SYNAPSES__
#define __H__UG__SYNAPSE_HANDLER__TYPES__SYNAPSES__

/* \defgroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace synapse_handler {
		/*!
		 * \brief synapses types
		 */
		enum SynapseType {
			EMPTY_SYNAPSE=0,
			ALPHA_SYNAPSE,
			EXP2_SYNAPSE,
			CUSTOM_SYNAPSE //!< could be used or another for custom synapses
		};
	} // namespace synapse_handler
} // namespace ug
//<! \}

#endif /// __H__UG__SYNAPSE_HANDLER__TYPES__SYNAPSES__
