/*!
 * \file utility.h
 * \brief utilities
 *
 *  Created on: Apr 10, 2015
 *      Author: stephan
 */

/// guards
#ifndef __H__UG__SYNAPSE_HANDLER__UTILITY__
#define __H__UG__SYNAPSE_HANDLER__UTILITY__

/// includes
#include "lib_disc/domain.h"

/*! \addtogroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug{
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
	} // end namespace sp
} // end namespace ug
///!< \}

#endif ///!<  __H__UG__SYNAPSE_HANDLER__UTILITY__
