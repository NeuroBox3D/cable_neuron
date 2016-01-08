/*!
 * \file utility.cpp
 *
 *  Created on: Apr 10, 2015
 *      Author: stephan
 */

/// includes
#include "utility.h"
#include "lib_disc/domain_traits.h"
#include "common/util/string_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/file_io/file_io_ugx.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"
#include "common/profiler/profiler.h"
#include "common/util/smart_pointer.h"

namespace ug{
namespace cable_neuron {
namespace synapse_handler {

	using namespace std;

	////////////////////////////////////////////////////////////////////////
	/// LoadDomainFromMultiGridInMemory
	////////////////////////////////////////////////////////////////////////
	template <typename TDomain>
		void LoadDomainFromGridInMemory(TDomain& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh) {
			LoadDomainFromGridInMemory(domain, grid, sh, 0);
	}


	////////////////////////////////////////////////////////////////////////
	/// LoadDomainFromMultiGridInMemory
	////////////////////////////////////////////////////////////////////////
	template <typename TDomain>
	void LoadDomainFromGridInMemory(TDomain& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int procId) {
		PROFILE_FUNC_GROUP("grid");
		domain.grid()->message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STARTS, procId));

		bool loadingMultiGrid = true;
		#ifdef UG_PARALLEL
			if((procId != -1) && (pcl::ProcRank() != procId))
				loadingMultiGrid = false;
		#endif

			if(loadingMultiGrid){
				if (!grid) {
					UG_THROW("ERROR in LoadDomainFromGridInMemory: Grid* was NULL.");
				} else {
					LoadGridFromMemory(domain, grid, sh);
				}
			}
			domain.grid()->message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STOPS, procId));
	}

	////////////////////////////////////////////////////////////////////////
	/// LoadGridFromMemory
	////////////////////////////////////////////////////////////////////////
	template <typename TDomain>
	void LoadGridFromMemory(TDomain& domain, const ug::MultiGrid* const gridIn, const ug::MGSubsetHandler* const shIn) {
		/// some checks
		if (!gridIn) {
			UG_THROW("ERROR in LoadGridFromMemory: MultiGrid* was NULL");
		} else if (!shIn) {
			UG_THROW("ERROR in LoadGridFromMemory: MGSubsetHandler* was NULL");
		}

		/// set grid, subset handler and assign grid
		ug::MultiGrid& grid = *domain.grid();
		grid = *gridIn;
		ug::MGSubsetHandler& mgsh = *domain.subset_handler();
		mgsh.assign_grid(grid);
		mgsh = *shIn;

		/// some debugging
		UG_LOGN("num_subsets (after): " << mgsh.num_subsets());
		UG_LOGN("num_levels (after): " << mgsh.num_levels());
	}

	/// explicit template instantiations
	template void LoadDomainFromGridInMemory<Domain1d>(Domain1d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh);
	template void LoadDomainFromGridInMemory<Domain2d>(Domain2d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh);
	template void LoadDomainFromGridInMemory<Domain3d>(Domain3d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh);
	template void LoadDomainFromGridInMemory<Domain1d>(Domain1d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int procId);
	template void LoadDomainFromGridInMemory<Domain2d>(Domain2d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int procId);
	template void LoadDomainFromGridInMemory<Domain3d>(Domain3d& domain, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int procId);

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

