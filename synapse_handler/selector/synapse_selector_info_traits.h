/*!
 * \file synapse_selector_info_traits.h
 *
 *  Created on: Apr 17, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__SYNAPSE_HANDLER__SYNAPSE_SELECTOR_INFO_TRAITS__
#define __H__UG__SYNAPSE_HANDLER__SYNAPSE_SELECTOR_INFO_TRAITS__

#include "../function/types.h"
/*! \addtogroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
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
} // namespace ug
//<! \}

#endif // __H__UG__SYNAPSE_HANDLER__SYNAPSE_SELECTOR_INFO_TRAITS__
