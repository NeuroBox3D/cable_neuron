/*
 * plugin_main.cpp
 *
 *  Created on: 03.11.2014
 *      Author: Lukas Reinhardt
 */

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"
#include "synapse_distributor.h"
#include <string>


using namespace std;
using namespace ug::bridge;
using namespace ug;


extern "C" UG_API void
InitUGPlugin_SynapseDistributor(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("SynapseDistributor/");
	string name = "SynapseDistributor";

	string constr_params("infile#outfile");

	typedef Attachment<vector<SynapseInfo> > AVSynapse;
	GlobalAttachments::declare_attachment<AVSynapse>("Synapses");

	typedef ug::SynapseDistributor TSD;

	reg->add_class_<TSD>(name, grp)
		.add_constructor<void (*)(string, string, bool)>("infile#outfile", "Initializes a SynapseDistributor", grp, "")
		.add_constructor<void (*)(SmartPtr<ug::Domain3d>, string, bool)>("dom#outfile", "Initializes a SynapseDistributor", grp, "")

		.add_method("remove_synapse", &TSD::remove_synapse, "e","","Removes synapse from edge",grp)
		.add_method("remove_all_synapses", &TSD::remove_all_synapses, "e","","Removes all synapses from edge",grp)
		.add_method("clear", static_cast<void (TSD::*)()>(&TSD::clear),"","","Removes all Synapses from grid",grp)
		.add_method("clear", static_cast<void (TSD::*)(int)>(&TSD::clear),"subsetIndex","","Removes all Synapses from subset",grp)
		.add_method("place_synapse", &TSD::place_synapse, "e", "", "Places a synapse.", grp)
		.add_method("place_synapses_uniform", static_cast<void (TSD::*)(size_t)>(&TSD::place_synapses_uniform), "", "", "Distributes synapses uniformly on grid.", grp)
		.add_method("place_synapses_uniform", static_cast<void (TSD::*)(int, size_t)>(&TSD::place_synapses_uniform), "", "", "Distributes synapses uniformly on subset.", grp)
		.add_method("place_synapses",static_cast<void (TSD::*)(std::vector<double>, size_t)>(&TSD::place_synapses), "p#subsetIndex", "", "Distributes synapses on all subsets s.t. a given density vector.", grp)
		.add_method("set_activation_timing",static_cast<void (TSD::*)(double, double, double, double)>(&TSD::set_activation_timing), "start_time#duration#start_time_dev#duration_dev", "", "Sets activity timing of distributed synapes..", grp)
		.add_method("degenerate_uniform",static_cast<void (TSD::*)(double)>(&TSD::degenerate_uniform),"","","Degenerates a certain percentage of synapses in the whole grid.",grp)
		.add_method("degenerate_uniform",static_cast<void (TSD::*)(double, int)>(&TSD::degenerate_uniform),"","","Degenerates a certain percentage of synapses in the given subset.",grp)

		.add_method("num_synapses", static_cast<size_t (TSD::*)()>(&TSD::num_synapses), "", "", "Returns global number of synapses", grp)
		.add_method("num_synapses", static_cast<size_t (TSD::*)(int)>(&TSD::num_synapses), "", "", "Returns number of synapses in specified subset", grp)
		.add_method("num_active_synapses", static_cast<size_t (TSD::*)(number)>(&TSD::num_active_synapses), "", "", "Returns global number of synapses at the specific time", grp)
		.add_method("num_active_synapses", static_cast<size_t (TSD::*)(number, int)>(&TSD::num_active_synapses), "", "", "Returns number of synapses at the specific in specified subset", grp)

		.add_method("activity_info",&TSD::activity_info, "", "","Prints start and end time for each synapse", grp)

		.add_method("print_status",&TSD::print_status, "t", "","prints synapse status of grid", grp)
		.add_method("get_last_message",&TSD::get_last_message,"","string","Returns last Message",grp)
		.add_method("export_grid", &TSD::export_grid,"","bool","Saves changes to disk",grp)
		.add_method("get_grid",&TSD::get_grid,"","Grid*","Pointer to current grid object",grp)
		.add_method("get_subset_handler",&TSD::get_subset_handler,"","SubsetHandler*","Pointer to current subsethandler object",grp)

		.set_construct_as_smart_pointer(true)
		;

}
