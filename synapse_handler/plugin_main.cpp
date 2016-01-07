/*!
 * \file plugin_main.cpp
 * \brief register synapse_handler plugin
 *
 *  Created on: Jan 19, 2015
 *      Author: stephan
 */

// includes
#include <string>
#include <map>
#include <vector>

#include <bridge/util.h>
#include <bridge/util_domain_dependent.h>
#include <common/error.h>
#include <common/ug_config.h>
#include <common/log.h>
#include <registry/registry.h>
#include <registry/error.h>
#include <lib_disc/domain.h>
#include <lib_grid/lib_grid.h>

#include "synapse_handler.h"
#include "util/utility.h"
#include "function/types.h"
#include "grid/synapse_info_io_traits.h"
#include "grid/synapse_info.h"

namespace ug {
namespace synapse_handler {

using namespace ug::bridge;
using namespace std;

/*!
 * \defgroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 *
 */

/// the functionality which is to be registered
struct Functionality
{
	/*!
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * \param reg registry
	 * \param parentGroup group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain(Registry& reg, string grp) {
		string suffix = GetDomainSuffix<TDomain>();
		string tag = GetDomainTag<TDomain>();

		/// implementations Vec1d
		typedef ISynapseHandler<TDomain> TISH;
		std::string name = std::string("ISynapseHandler").append(suffix);
		reg.add_class_<TISH>(name, grp)
			.add_method("name", &TISH::name);
		reg.add_class_to_group(name, "ISynapseHandler", tag);

#ifdef UG_FOR_LUA
		typedef NETISynapseHandler<TDomain> TNETISH;
		name = std::string("NETISynapseHandler").append(suffix);
		reg.add_class_<TNETISH, TISH>(name, grp)
			.template add_constructor<void (*)()> ()
			.add_method("set_vmdisc", &TNETISH::set_vmdisc)
			.add_method("set_presyn_subset", &TNETISH::set_presyn_subset)
			.add_method("set_activation_timing",
				static_cast<void (TNETISH::*)(number, number, number, number)>(&TNETISH::set_activation_timing),
				"", "start_time#duration#start time deviation#duration deviation", "")
			.add_method("set_activation_timing",
				static_cast<void (TNETISH::*)(number, number, number, number, number)>(&TNETISH::set_activation_timing),
				"", "start_time#duration#start time deviation#duration deviation#peak conductivity", "")
			.add_method("set_activation_timing",
				static_cast<void (TNETISH::*)(number, number, number, number, number, bool)>(&TNETISH::set_activation_timing),
				"", "start_time#duration#start time deviation#duration deviation#peak conductivity#seed", "")
			//.add_method("set_custom_diameter", &TNETISH::set_custom_diameter, "", "subsets#diameter", "")
			.add_method("update_presyn", &TNETISH::update_presyn)
			.add_method("print_synapse_statistics", &TNETISH::print_synapse_statistics, "", "soma subset index", "")
			.add_method("write_activity_to_file", &TNETISH::write_activity_to_file, "", "file base name#time", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NETISynapseHandler", tag);
#endif

#ifdef SH_SYNAPSE_DISTRIBUTOR_ENABLED
		typedef SynapseDistributorSynapseHandler<TDomain> TSDSH;
		name = std::string("SDSynapseHandler").append(suffix);
		reg.add_class_<TSDSH, TISH>(name, grp)
			.add_method("set_sd", &TSDSH::set_sd)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SDSynapseHandler", tag);
#endif

		/// utility
		reg.add_function("LoadDomainFromGridInMemory", static_cast<void (*)(TDomain&, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int)>(
							 &LoadDomainFromGridInMemory<TDomain>));
		reg.add_function("LoadDomainFromGridInMemory", static_cast<void (*)(TDomain&, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh)>(
							 &LoadDomainFromGridInMemory<TDomain>));
	}

	/*!
	 * Function called for the registration of Domain and Algebra independent parts.
	 * All Functions and Classes not depending on Domain and Algebra
	 * are to be placed here when registering.
	 *
	 * \param reg		registry
	 * \param grp		group for sorting of functionality
	 */
	static void Common(Registry& reg, string grp)
	{
	}

// end of functionality which is to be exported
};
/// \}


// end namespace sh
}


/// \addtogroup sp_plugin
extern "C" void
InitUGPlugin_SynapseHandler(bridge::Registry* reg, std::string grp)
{
	typedef Attachment<uint> APresynInd;
	typedef Attachment<std::vector<SynapseInfo> >AVSynapse;
	GlobalAttachments::declare_attachment<APresynInd>("presyn_index", true);
	GlobalAttachments::declare_attachment<AVSynapse>("Synapses", false);

	grp.append("/SynapseHandler");
	typedef synapse_handler::Functionality SynapseHandlerFunctionality;
	try
	{
		//RegisterCommon<SynapseHandlerFunctionality>(*reg,grp);
		RegisterDomainDependent<SynapseHandlerFunctionality>(*reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

// end namespace ug
}
