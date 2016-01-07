/*!
 * \file synapse_factory.h
 * \brief creates a synapse
 *
 *  Created on: Apr 21, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__SYNAPSE_HANDLER__SYNAPSE_FACTORY__
#define __H__UG__SYNAPSE_HANDLER__SYNAPSE_FACTORY__

/// includes
#include "synapse.h"

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
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
} // namespace ug
//<! \}

#endif // __H__UG__SYNAPSE_HANDLER__SYNAPSE_FACTORY__
