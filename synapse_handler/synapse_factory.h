/*!
 * \file synapse_factory.h
 * \brief creates a synapse
 *
 *  Created on: Apr 21, 2015
 *      Author: stephan
 */

/// guard
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_FACTORY_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_FACTORY_H__

/// includes
#include "synapse.h"

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {
	/*!
	 * \brief synapse factory interface
	 */
	class ISynapseFactory {
	public:
		/*!
		 * \brief creates a synapse
		 */
		virtual ISynapse* create_synapse() = 0;

		/*!
		 * \brief mandatory vdtor
		 */
		virtual ~ISynapseFactory() {
		}
	};
} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug
//<! \}

#endif // __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__SYNAPSE_FACTORY_H__
