/*!
 * \file utility.h
 * \brief utilities
 *
 *  Created on: Apr 10, 2015
 *      Author: stephan
 */

/// guards
#ifndef __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__UTIL__UTILITY_H__
#define __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__UTIL__UTILITY_H__

/// includes
#include "lib_disc/domain.h"

/*! \addtogroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace cable_neuron {
namespace synapse_handler {
	/*!
	 * \brief Loads a domain from an existing grid in memory
	 * \tparam[out] domain
	 * \param[in] grid
	 */
	template <typename TDomain>
	void LoadDomainFromGridInMemory(TDomain& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh);

	/*!
	 * \brief Loads a domain from an ug::MultiGrid in memory with given procId
	 * \tparam[out] domain
	 * \param[in] grid
	 * \param[in] procId
	 */
	template <typename TDomain>
	void LoadDomainFromGridInMemory(TDomain& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int procId);

	/*!
	 * \brief set pointer to grid
	 * \tparam[out] domain
	 * \param[in] gridIn
	 */
	template <typename TDomain>
	void LoadGridFromMemory(TDomain& domain, const ug::MultiGrid* const gridIn, const ug::MGSubsetHandler* const sh);
} // namespace synapse_handler
} // namespace cable_neuron
} // end namespace ug
///!< \}

#endif ///!<  __UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLER__UTIL__UTILITY_H__
