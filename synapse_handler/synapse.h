/*!
 * \file synapse.h
 * \brief synapse interface
 *
 *  Created on: Mar 23, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__SYNAPSE_HANDLER__SYNAPSE__
#define __H__UG__SYNAPSE_HANDLER__SYNAPSE__

/* \defgroup sh_plugin Synapse Handler plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace synapse_handler {
		/*!
		 * \brief synapse interface
		 */
		class ISynapse {
		public:
			/*!
			 * \brief return the current
			 * \param[in] vm membrane potential
			 * \param[in] time current time
			 */
			virtual number current(const number& t) const = 0;

			/*!
			 * \brief mandatory virtual dtor
			 */
			virtual ~ISynapse() {

			}
		};
	} // namespace synapse_handler
} // namespace ug
//<! \}

#endif // __H__UG__SYNAPSE_HANDLER__SYNAPSE__
