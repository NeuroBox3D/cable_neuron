/*
 * plugin_main.cpp
 *
 *  Created on: 2014-09-02
 *  Author: mbreit, pgottmann
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "cable_disc/cable_equation.h"
#include "cable_disc/ElemDiscHH_base.h"
#include "cable_disc/ElemDiscHH_fv1.h"
#include "cable_disc/ElemDiscHH_Nernst_fv1.h"
#include "cable_disc/ElemDiscHH_Nernst_neuron_fv1.h"

#include "lib_grid/global_attachments.h" // global attachments
#include "lib_disc/function_spaces/grid_function.h"

// cable equation includes
#include "membrane_transport/cable_membrane_transport_interface.h"

// solver includes
#include "util/cable_ass_tuner.h"
#include "util/order.h"

// implemented membrane transporters
#include "membrane_transport/channel_hh.h"
#include "membrane_transport/leakage.h"
#include "membrane_transport/vdcc_bg.h"
#include "membrane_transport/ion_leakage.h"
#include "membrane_transport/na_k_pump.h"
#include "membrane_transport/ncx.h"
#include "membrane_transport/pmca.h"

// converted membrane transporters
#ifdef CONVERTED_TRANSPORT_ENABLED
	#include "membrane_transport/nmodl_converter/converted/includefile.h"
#endif

// synapse handler
#include "synapse_handler/synapse_handler.h"
#include "synapse_handler/util/utility.h"
#include "synapse_handler/function/types.h"
#include "synapse_handler/grid/synapse_info_io_traits.h"
#include "synapse_handler/grid/synapse_info.h"

//split_synapse handler
//#include "split_synapse_handler/IPreSynapse.h"

// synapse distributor
#include "synapse_distributor/synapse_distributor.h"

// utility
#include "util/functions.h"


using namespace std;
using namespace ug::bridge;

namespace ug {
namespace cable_neuron {



/**
 *  \defgroup plugin_cable_neuron Plugin cable_neuon
 *  \ingroup plugins_experimental
 *  This is a plugin for cable equation functionality.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

	/**
	 * Function called for the registration of Domain dependent parts
	 * of the plugin. All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain(Registry& reg, string grp)
	{
		static const int dim = TDomain::dim;
		string suffix = GetDomainSuffix<TDomain>();
		string tag = GetDomainTag<TDomain>();


		// /////////////////////////////////////////
		// //////// synapse handler ////////////////
		// /////////////////////////////////////////

		{
			using namespace synapse_handler;

			// synapse handler interface
			typedef ISynapseHandler<TDomain> TISH;
			std::string name = std::string("ISynapseHandler").append(suffix);
			reg.add_class_<TISH>(name, grp)
				.add_method("name", &TISH::name);
			reg.add_class_to_group(name, "ISynapseHandler", tag);

			// NETI synapse handler
			typedef NETISynapseHandler<TDomain> TNETISH;
			name = std::string("NETISynapseHandler").append(suffix);
			reg.add_class_<TNETISH, TISH>(name, grp)
				.template add_constructor<void (*)()> ()
				.add_method("set_ce_object", &TNETISH::set_ce_object)
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
				//.add_method("update_presyn", &TNETISH::update_presyn)	// no, handled internally
				.add_method("print_synapse_statistics", &TNETISH::print_synapse_statistics, "", "soma subset index", "")
				.add_method("write_activity_to_file", &TNETISH::write_activity_to_file, "", "file base name#time", "")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NETISynapseHandler", tag);

			// SynapseDistributor synapse handler
			typedef SynapseDistributorSynapseHandler<TDomain> TSDSH;
			name = std::string("SDSynapseHandler").append(suffix);
			reg.add_class_<TSDSH, TISH>(name, grp)
				.add_method("set_sd", &TSDSH::set_sd)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "SDSynapseHandler", tag);

			// utility
			reg.add_function("LoadDomainFromGridInMemory",
				static_cast<void (*)(TDomain&, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh, int)>
					(&LoadDomainFromGridInMemory<TDomain>));
			reg.add_function("LoadDomainFromGridInMemory",
				static_cast<void (*)(TDomain&, const ug::MultiGrid* const grid, const ug::MGSubsetHandler* const sh)>
					(&LoadDomainFromGridInMemory<TDomain>));
		}



		// ////////////////////////////////////
		// //////// cable disc ////////////////
		// ////////////////////////////////////

		// Kabel Diff Base
		{
			typedef ElemDiscHH_Base<TDomain> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("ElemDiscHH_Base").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
					//.add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)");
						.add_method("set_spec_capa", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_spec_capa), "", "Spec Capa")
						.add_method("set_spec_capa", static_cast<void (T::*)(number)>(&T::set_spec_capa), "", "Spec Capa")
			#ifdef UG_FOR_LUA
						.add_method("set_spec_capa", static_cast<void (T::*)(const char*)>(&T::set_spec_capa), "", "Spec Capa");
			#else
			;
			#endif
			reg.add_class_to_group(name, "ElemDiscHH_Base", tag);
		}



		// Kabel Diff FV1 with constant reversal potential of K and Na
		{
			typedef ElemDiscHH_FV1<TDomain> T;
			typedef ElemDiscHH_Base<TDomain> TBase;
			string name = string("ElemDiscHH_FV1").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*,const char*)>("Function(s)#Subset(s)")
				.add_method("set_injection", &T::set_injection)
				.add_method("set_diameter", &T::set_diameter)
				.add_method("set_spec_res", &T::set_spec_res)
				.add_method("set_consts", &T::set_consts)
				.add_method("set_rev_pot", &T::set_rev_pot)
				.add_method("set_accuracy", &T::set_accuracy)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ElemDiscHH_FV1", tag);
		}



		// Kabel Diff FV1 with dynamic calculated reversal potential of K and Na (Nernst-Equatation)
		{
			typedef ElemDiscHH_Nernst_FV1<TDomain> T;
			typedef ElemDiscHH_Base<TDomain> TBase;
			string name = string("ElemDiscHH_Nernst_FV1").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*,const char*)>("Function(s)#Subset(s)")
				.add_method("set_injection", &T::set_injection)
				.add_method("set_diameter", &T::set_diameter)
				.add_method("set_spec_res", &T::set_spec_res)
				.add_method("set_consts", &T::set_consts)
				.add_method("set_nernst_consts", &T::set_nernst_consts)
				.add_method("set_accuracy", &T::set_accuracy)
				.add_method("set_diff_Na", &T::set_diff_Na)
				.add_method("set_diff_K", &T::set_diff_K)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ElemDiscHH_Nernst_FV1", tag);
		}

		// Kabel Diff FV1 with dynamic calculated reversal potential of K and Na (Nernst-Equatation)
		// Gating functions from neuron
		{
			typedef ElemDiscHH_Nernst_neuron_FV1<TDomain> T;
			typedef ElemDiscHH_Base<TDomain> TBase;
			string name = string("ElemDiscHH_Nernst_neuron_FV1").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*,const char*)>("Function(s)#Subset(s)")
				.add_method("set_injection", &T::set_injection)
				.add_method("set_diameter", &T::set_diameter)
				.add_method("set_spec_res", &T::set_spec_res)
				.add_method("set_consts", &T::set_consts)
				.add_method("set_nernst_consts", &T::set_nernst_consts)
				.add_method("set_accuracy", &T::set_accuracy)
				.add_method("set_diff_Na", &T::set_diff_Na)
				.add_method("set_diff_K", &T::set_diff_K)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ElemDiscHH_Nernst_neuron_FV1", tag);
		}



		// Channel Interface Base (virtual class)
		{
			typedef ICableMembraneTransport<TDomain> T;
			string name = string("ICableMembraneTransport").append(suffix);
			reg.add_class_<T>(name, grp);
				//.add_method("init", &T::init)
				//.add_method("update_gating", &T::update_gating);
			reg.add_class_to_group(name, "ICableMembraneTransport", tag);
		}


		// HH
		{
			typedef ChannelHH<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("ChannelHH").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number)>(&T::set_conductances), "",
						"K conductance (S/m^2)| default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number, const char*)>(&T::set_conductances), "",
						"K conductance (S/m^2) | default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3#"
						"subset(s) as C-type string",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number, const std::vector<std::string>&)>(&T::set_conductances), "",
						"K conductance (S/m^2) | default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3#"
						"subset(s) as vector of string",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_log_mGate", &T::set_log_mGate)
				.add_method("set_log_nGate", &T::set_log_nGate)
				.add_method("set_log_hGate", &T::set_log_hGate)
				//.add_method("current", /*static_cast<void (TBase::*) (Vertex*, std::vector<>&)> (&T::current) /*, "","", "doing flux")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelHH", tag);
		}


		// HH-with-Nernst
		{
			typedef ChannelHHNernst<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("ChannelHHNernst").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number)>(&T::set_conductances), "",
						"K conductance (S/m^2) | default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number, const char*)>(&T::set_conductances), "",
						"K conductance (S/m^2) | default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3#"
						"subset(s) as C-type string",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_conductances",  static_cast<void (T::*)(number, number, const std::vector<std::string>&)>(&T::set_conductances), "",
						"K conductance (S/m^2) | default | value=3.6e2#"
						"Na conductance (S/m^2) | default | value=1.2e3#"
						"subset(s) as vector of string",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_log_mGate", &T::set_log_mGate)
				.add_method("set_log_nGate", &T::set_log_nGate)
				.add_method("set_log_hGate", &T::set_log_hGate)
				//.add_method("current", /*static_cast<void (TBase::*) (Vertex*, std::vector<>&)> (*/&T::current) /*, "","", "doing flux")*/
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelHHNernst", tag);
		}


		// charge leakage
		{
			typedef ChannelLeak<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("ChannelLeak").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_cond", static_cast<void (T::*)(number)>(&T::set_cond),
							"", "leak conductance (S/m^2) | default | value=1.0",
							"sets leak conductance for leak channel")
				.add_method("set_cond", static_cast<void (T::*)(number, const char*)>(&T::set_cond),
							"", "leak conductance (S/m^2) | default | value=1.0 # subset(s) as C-type string",
							"sets leak conductance for leak channel")
				.add_method("set_cond", static_cast<void (T::*)(number, const std::vector<std::string>&)>(&T::set_cond),
							"", "leak conductance (S/m^2) | default | value=1.0 # subset(s) as vector of strings",
							"sets leak conductance for leak channel")
				.add_method("set_rev_pot",  static_cast<void (T::*)(number)>(&T::set_rev_pot), "",
							"leakage equilibrium potential (V) | default | value=-0.065",
							"sets leakage equilibrium potential")
				.add_method("set_rev_pot",  static_cast<void (T::*)(number, const char*)>(&T::set_rev_pot),
							"", "leakage equilibrium potential (V) | default | value=-0.065 # subset(s) as C-type string",
							"sets leakage equilibrium potential")
				.add_method("set_rev_pot",  static_cast<void (T::*)(number, const std::vector<std::string>&)>(&T::set_rev_pot),
							"", "leakage equilibrium potential (V) | default | value=-0.065 # subset(s) as vector of strings",
							"sets leakage equilibrium potential")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelLeak", tag);
		}


		// ion leakage
		{
			typedef IonLeakage<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("IonLeakage").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_perm", &T::set_perm, "", "flux at rest (mol/(m^2*s))#inner concentration at rest (mM)# outer concentration at rest (mM)#potential at rest (V)",
						    "tunes the channel to equilibrate fluxes at resting conditions")
				.add_method("set_leaking_quantity", &T::set_leaking_quantity, "", "", "")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "IonLeakage", tag);
		}


		// VDCC BG
		{
			typedef VDCC_BG_cable<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("VDCC_BG_cable").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_log_mGate" , &T::set_log_mGate)
				.add_method("set_log_hGate" , &T::set_log_hGate)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VDCC_BG_cable", tag);
		}

		// PMCA
		{
			typedef PMCA_cable<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("PMCA_cable").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_max_flux", &T::set_max_flux, "", "max. flux density (mol/(m^2*s)) | default | value=8.5e-9", "sets maximal flux density")
				.add_method("set_kd", &T::set_kd, "", "K_D value inner [Ca] (mM) | default | value=6.0e-5", "sets K_D value for inner [Ca]")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "PMCA_cable", tag);
		}

		// NCX
		{
			typedef NCX_cable<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("NCX_cable").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_max_flux", &T::set_max_flux, "", "max. flux density (mol/(m^2*s)) | default | value=3.75e-8", "sets maximal flux density")
				.add_method("set_kd", &T::set_kd, "", "K_D value inner [Ca] (mM) | default | value=1.8e-3", "sets K_D value for inner [Ca]")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NCX_cable", tag);
		}


		// Na/K pump
		{
			typedef Na_K_Pump<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("Na_K_Pump").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_max_flux", &T::set_max_flux, "", "max. flux density (mol/(m^2*s)) | default | value=3.6e-2", "sets maximal flux density")
				.add_method("set_K_K", &T::set_K_K, "", "K_D value outer [K] (mM) | default | value=1.37", "sets K_D value for outer [K]")
				.add_method("set_K_Na", &T::set_K_Na, "", "K_D value inner [Na] (mM) | default | value=5.74", "sets K_D value for inner [Na]")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Na_K_Pump", tag);
		}

		// CableEquation discretization class
		{
			typedef CableEquation<TDomain> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("CableEquation").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*)>
					("Subset(s)")
				.template add_constructor<void (*)(const char*, bool)>
					("Subset(s)#with ion concentrations?")
				.template add_constructor<void (*)(const char*, bool, number)>
					("Subset(s)#with ion concentrations?#InitTime")
				.add_method("set_diameter", static_cast<void (T::*)(number)>(&T::set_diameter),
						"", "diameter (m) | default | value=1e-6", "sets a new diameter")

				.add_method("set_spec_res", static_cast<void (T::*)(number)>(&T::set_spec_res),
						"", "specific resistance (Ohm m) | default | value=1.0", "sets a new specific resistance")
				.add_method("set_spec_cap", static_cast<void (T::*)(number)>(&T::set_spec_cap),
						"", "specific capacitance (F/m^2) | default | value=1e-2", "sets a new specific capacitance")

				.add_method("set_k_out", static_cast<void (T::*)(number)>(&T::set_k_out),
						"", "extracellular [K] (mM) | default | value=4.0", "sets extracellular [K]")
				.add_method("set_na_out", static_cast<void (T::*)(number)>(&T::set_na_out),
						"", "extracellular [Na] (mM) | default | value=150.0", "sets extracellular [Na]")
				.add_method("set_ca_out", static_cast<void (T::*)(number)>(&T::set_ca_out),
						"", "extracellular [Ca] (mM) | default | value=1.5", "sets extracellular [Ca]")

				.add_method("set_rev_pot_na", static_cast<void (T::*)(number)>(&T::set_rev_pot_na),
						"", "Na reversal potential (V) | default | value=0.06", "sets reversal potential for Na")
				.add_method("set_rev_pot_k", static_cast<void (T::*)(number)>(&T::set_rev_pot_k),
						"", "K reversal potential (V) | default | value=-0.09", "sets reversal potential for K")
				.add_method("set_rev_pot_ca", static_cast<void (T::*)(number)>(&T::set_rev_pot_ca),
						"", "Ca reversal potential (V) | default | value=0.14", "sets reversal potential for Ca")

				.add_method("set_temperature", static_cast<void (T::*)(number)>(&T::set_temperature),
						"", "temperature (K) | default | value=310", "sets new temperature")
				.add_method("set_temperature_celsius", static_cast<void (T::*)(number)>(&T::set_temperature_celsius),
						"", "temperature (°C) | default | value=37", "sets new temperature")

				.add_method("set_diff_coeffs", static_cast<void (T::*)(const std::vector<number>&)> (&T::set_diff_coeffs), "",
						"diffusion coefficients of K, Na and Ca (m^2/s)", "sets diffusion coefficients")

				.add_method("add", &T::add)
				.add_method("set_influx", static_cast<void (T::*)(number, number, number, number, number, number)>(&T::set_influx), "",
						"current (A) | default | value=1e-9 #"
						"x-coordinate of influx position (m) | default | 0.0 #"
						"y-coordinate of influx position (m) | default | 0.0 #"
						"z-coordinate of influx position (m) | default | 0.0 #"
						"start time (s) | default | 0 #"
						"duration (s) | default | 0 ",
						"sets position, duration and current strength of an influx")
				.add_method("write_states_for_position", &T::write_states_for_position)
				.add_method("set_output_point_and_path", &T::set_output_point_and_path)
				.add_method("set_influx_subset", &T::set_influx_subset)
#ifdef UG_CPU_1
				.add_method("estimate_cfl_cond", &T::template estimate_cfl_cond<CPUAlgebra::vector_type>)
#endif
				.add_method("set_synapse_handler", &T::set_synapse_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "CableEquation", tag);
		}

#ifdef CONVERTED_TRANSPORT_ENABLED
		#include "membrane_transport/nmodl_converter/converted/channels.cpp"
#endif

		// Cuthill-McKee ordering
		{
			reg.add_function("order_cuthillmckee", &order_cuthillmckee<TDomain>, grp.c_str(),
							 "", "approxSpace", "vertex ordering for solver optimization");
		}

		// checks for acyclicity, presynaptic indices, both
		{
			reg.add_function("is_acyclic", static_cast<bool (*) (SmartPtr<TDomain>)>(&is_acyclic<TDomain>), grp.c_str(),
							 "", "domain", "Checks whether given domain is acyclic.");
			reg.add_function("is_acyclic", static_cast<bool (*) (SmartPtr<TDomain>, int)>(&is_acyclic<TDomain>), grp.c_str(),
							 "", "domain, verbosity", "Checks whether given domain is acyclic.");
			reg.add_function("check_presyn_indices", static_cast<int (*) (SmartPtr<TDomain>)>(&check_presyn_indices<TDomain>), grp.c_str(),
							 "", "domain", "Checks whether presynaptic indices are in order.");
			reg.add_function("check_presyn_indices", static_cast<int (*) (SmartPtr<TDomain>, int)>(&check_presyn_indices<TDomain>), grp.c_str(),
							 "", "domain, verbosity", "Checks whether presynaptic indices are in order.");
			reg.add_function("check_domain", static_cast<int (*) (SmartPtr<TDomain>)>(&check_domain<TDomain>), grp.c_str(),
							 "", "domain", "Checks whether given domain is sound.");
			reg.add_function("check_domain", static_cast<int (*) (SmartPtr<TDomain>, int)>(&check_domain<TDomain>), grp.c_str(),
							 "", "domain, verbosity", "Checks whether given domain is sound.");
		}
	}

	/**
	 * Function called for the registration of Domain and Algebra dependent parts.
	 * All Functions and Classes depending on both Domain and Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Domain and Algebra types, based on the current build options.
	 *
	 * @param reg		registry
	 * @param grp		group for sorting of functionality
	 */
	template <typename TDomain, typename TAlgebra>
	static void DomainAlgebra(bridge::Registry& reg, string grp)
	{
		string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
		string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

		// CableAssTuner
		{
			typedef CableAssTuner<TDomain, TAlgebra> T;
			string name = string("CableAssTuner");
			reg.add_class_<T>(name+suffix, grp)
				.template add_constructor<void (*)(SmartPtr<DomainDiscretization<TDomain, TAlgebra> >,
												   SmartPtr<ApproximationSpace<TDomain> >)>("domain disc")
				.add_method("remove_ghosts_from_assembling_iterator", &T::template remove_ghosts_from_assembling_iterator<1>)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name+suffix, name, tag);
		}
	}

	/**
	 * Function called for the registration of Dimension dependent parts.
	 * All Functions and Classes depending on the Dimension
	 * are to be placed here when registering. The method is called for all
	 * available Dimension types, based on the current build options.
	 *
	 * @param reg		registry
	 * @param grp		group for sorting of functionality
	 */
	template <int dim>
	static void Dimension(Registry& reg, string grp)
	{
		string suffix = GetDimensionSuffix<dim>();
		string tag = GetDimensionTag<dim>();
	}

	/**
	 * Function called for the registration of Algebra dependent parts.
	 * All Functions and Classes depending on Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Algebra types, based on the current build options.
	 *
	 * @param reg		registry
	 * @param grp		group for sorting of functionality
	 */
	template <typename TAlgebra>
	static void Algebra(Registry& reg, string grp)
	{
		string suffix = GetAlgebraSuffix<TAlgebra>();
		string tag = GetAlgebraTag<TAlgebra>();
	}

	/**
	 * Function called for the registration of Domain and Algebra independent parts.
	 * All Functions and Classes not depending on Domain and Algebra
	 * are to be placed here when registering.
	 *
	 * @param reg		registry
	 * @param grp		group for sorting of functionality
	 */
	static void Common(Registry& reg, string grp)
	{
		// /////////////////////////////////////////////
		// //////// synapse distributor ////////////////
		// /////////////////////////////////////////////

		// TODO: might better be registered (and implemented) in a domain-dependent manner
		string name = "SynapseDistributor";

		string constr_params("infile#outfile");


		typedef SynapseDistributor TSD;

		reg.add_class_<TSD>(name, grp)
			.add_constructor<void (*)(string, string, bool)>("infile#outfile", "Initializes a SynapseDistributor", grp, "")
			.add_constructor<void (*)(SmartPtr<ug::Domain3d>, string, bool)>("dom#outfile", "Initializes a SynapseDistributor", grp, "")

			.add_method("remove_synapse", &TSD::remove_synapse, "e","","Removes synapse from edge",grp)
			.add_method("remove_all_synapses", &TSD::remove_all_synapses, "e","","Removes all synapses from edge",grp)
			.add_method("clear", static_cast<void (TSD::*)()>(&TSD::clear),"","","Removes all Synapses from grid",grp)
			.add_method("clear", static_cast<void (TSD::*)(int)>(&TSD::clear),"subsetIndex","","Removes all Synapses from subset",grp)
			.add_method("place_synapse", &TSD::place_synapse, "e", "", "Places a synapse.", grp)
			.add_method("place_synapses_uniform", static_cast<void (TSD::*)(size_t)>(&TSD::place_synapses_uniform), "", "", "Distributes synapses uniformly on grid.", grp)
			.add_method("place_synapses_uniform", static_cast<void (TSD::*)(int, size_t)>(&TSD::place_synapses_uniform), "", "", "Distributes synapses uniformly on subset.", grp)
			.add_method("place_synapses_uniform", static_cast<void (TSD::*)(const char*, number)>(&TSD::place_synapses_uniform), "", "", "Distributes synapses uniformly on subset s.t. given density of synapses.", grp)
			.add_method("place_synapses",static_cast<void (TSD::*)(std::vector<number>, size_t)>(&TSD::place_synapses), "p#subsetIndex", "", "Distributes synapses on all subsets s.t. a given density vector.", grp)
			.add_method("set_activation_timing",static_cast<void (TSD::*)(number, number, number, number)>(&TSD::set_activation_timing), "start_time#duration#start_time_dev#duration_dev", "", "Sets activity timing of distributed synapes..", grp)
			.add_method("degenerate_uniform",static_cast<void (TSD::*)(number)>(&TSD::degenerate_uniform),"","","Degenerates a certain percentage of synapses in the whole grid.",grp)
			.add_method("degenerate_uniform",static_cast<void (TSD::*)(number, int)>(&TSD::degenerate_uniform),"","","Degenerates a certain percentage of synapses in the given subset.",grp)
			.add_method("degenerate_uniform",static_cast<void (TSD::*)(number, const char*)>(&TSD::degenerate_uniform),"","","Degenerates a certain percentage of synapses in the given subset.",grp)

			.add_method("get_subset_length",static_cast<number (TSD::*)(int)>(&TSD::get_subset_length),"","number","Calculate and return length of specified subset in micrometer",grp)
			.add_method("get_subset_length",static_cast<number (TSD::*)(const char*)>(&TSD::get_subset_length),"","number","Calculate and return length of specified subset in micrometer",grp)
			.add_method("num_synapses", static_cast<size_t (TSD::*)()>(&TSD::num_synapses), "", "", "Returns global number of synapses", grp)
			.add_method("num_synapses", static_cast<size_t (TSD::*)(int)>(&TSD::num_synapses), "", "", "Returns number of synapses in specified subset", grp)
			.add_method("num_synapses", static_cast<size_t (TSD::*)(const char*)>(&TSD::num_synapses), "", "", "Returns number of synapses in specified subset", grp)
			.add_method("num_active_synapses", static_cast<size_t (TSD::*)(number)>(&TSD::num_active_synapses), "", "", "Returns global number of synapses at the specific time", grp)
			.add_method("num_active_synapses", static_cast<size_t (TSD::*)(number, int)>(&TSD::num_active_synapses), "", "", "Returns number of synapses at the specific in specified subset", grp)

			.add_method("activity_info",&TSD::activity_info, "", "","Prints start and end time for each synapse", grp)

			.add_method("print_status",&TSD::print_status, "t", "","prints synapse status of grid", grp)
			.add_method("get_last_message",&TSD::get_last_message,"","string","Returns last Message",grp)
			.add_method("set_outfile",&TSD::set_outfile,"","void","Sets output filename",grp)
			.add_method("export_grid",static_cast<bool (TSD::*)()>(&TSD::export_grid),"","bool","Saves changes to disk",grp)
			.add_method("export_grid",static_cast<bool (TSD::*)(string)>(&TSD::export_grid),"","bool","Saves changes to disk",grp)
			.add_method("get_grid",&TSD::get_grid,"","Grid*","Pointer to current grid object",grp)
			.add_method("get_subset_handler",&TSD::get_subset_handler,"","SubsetHandler*","Pointer to current subsethandler object",grp)

			.set_construct_as_smart_pointer(true);

	}

}; // end Functionality


// end group plugin_cable_neuron
/// \}

} // end namespace cable_neuron



/// This function is called when the plugin is loaded.
extern "C" void
InitUGPlugin_cable_neuron(Registry* reg, string grp)
{
	grp.append("/cable_neuron");

	// declare global grid attachments for diameter, synapses and presynaptic indices
	typedef ANumber ADiameter;
	typedef Attachment<uint> APresynInd;
	typedef Attachment<std::vector<SynapseInfo> >AVSynapse;

	GlobalAttachments::declare_attachment<ADiameter>("diameter", true);
	GlobalAttachments::declare_attachment<APresynInd>("presyn_index", true);
	GlobalAttachments::declare_attachment<AVSynapse>("Synapses", false);

	typedef cable_neuron::Functionality Functionality;

	try
	{
		RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}	// namespace ug

