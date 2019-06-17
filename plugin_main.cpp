/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Pascal Gottmann
 * Creation date: 2014-09-02
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "cable_disc/cable_equation.h"
#include "cable_disc/cable_equation_withOuterPot.h"
#include "cable_disc/implicit_active_cable_disc.h"
#include "cable_disc/implicit_active_cable_disc_base.h"
#include "cable_disc/implicit_active_cable_disc_nernst.h"

#include "lib_grid/global_attachments.h" // global attachments
#include "lib_disc/function_spaces/grid_function.h"
#ifdef UG_PARALLEL
#include "lib_disc/parallelization/domain_load_balancer.h"  // for DomainPartitioner
#endif

// configuration file for compile options
#include "cn_config.h"

// cable equation includes
#include "membrane_transport/cable_membrane_transport_interface.h"

// solver includes
#include "util/cable_ass_tuner.h"
#include "util/order.h"
#include "util/neuronal_topology_importer.h"

// implemented membrane transporters
#include "membrane_transport/channel_hh.h"
#include "membrane_transport/leakage.h"
#include "membrane_transport/vdcc_bg.h"
#include "membrane_transport/ion_leakage.h"
#include "membrane_transport/na_k_pump.h"
#include "membrane_transport/ncx.h"
#include "membrane_transport/pmca.h"
#include "membrane_transport/ka_golding01.h"
#include "membrane_transport/kdr_golding01.h"
#include "membrane_transport/nax_golding01.h"

// converted membrane transporters
#ifdef CONVERTED_TRANSPORT_ENABLED
	#include "membrane_transport/nmodl_converter/converted/includefile.h"
#endif

// converted bap modls
#ifdef CONVERTED_TRANSPORT_BAP_ENABLED
	#include "membrane_transport/nmodl_converter/converted/includefile_bAP.h"
#endif

// synapse handling
#include "synapse_handling/synapses/base_synapse.h"
#include "synapse_handling/synapse_info_io_traits.h"
#include "synapse_handling/synapses/onset_pre_synapse.h"
#include "synapse_handling/synapses/alpha_post_synapse.h"
#include "synapse_handling/synapses/threshold_pre_synapse.h"
#include "synapse_handling/synapses/exp2_post_synapse.h"
#include "synapse_handling/synapse_dealer.h"
#include "synapse_handling/synapse_distributor.h"
#include "synapse_handling/synapse_container.h"
#include "synapse_handling/synapse_handler.h"

// utility
#include "util/functions.h"

// partitioning
#ifdef UG_PARALLEL
#include "util/simple_cable_neuron_partitioner.h"
#endif
#ifdef NC_WITH_PARMETIS
#include "util/cable_neuron_unificator.h"
#endif

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


		// ////////////////////////////////////
		// //////// cable disc ////////////////
		// ////////////////////////////////////

		// fully implicit cable equation (base class)
		{
			typedef ImplicitActiveCableDiscBase<TDomain> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("ElemDiscHH_Base").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.add_method("set_spec_cap", static_cast<void (T::*)(number)>(&T::set_spec_cap), "", "specific capacitance (in F/m^2)")
				.add_method("set_spec_cap", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_spec_cap), "", "specific capacitance (in F/m^2) function")
			#ifdef UG_FOR_LUA
				.add_method("set_spec_cap", static_cast<void (T::*)(const char*)>(&T::set_spec_cap), "", "specific capacitance (in F/m^2) function name")
			#endif

				.add_method("set_spec_res", static_cast<void (T::*)(number)>(&T::set_spec_res), "", "specific resistance (in Ohm m)")
				.add_method("set_spec_res", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_spec_res), "", "specific resistance (in Ohm m) function")
			#ifdef UG_FOR_LUA
				.add_method("set_spec_res", static_cast<void (T::*)(const char*)>(&T::set_spec_res), "", "specific resistance (in Ohm m) function name")
			#endif

				.add_method("set_conductances", static_cast<void (T::*)(number, number, number)>(&T::set_conductances), "", "conductances for K+, Na+ and leakage (in S/m^2)")
				.add_method("set_conductances", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, SmartPtr<CplUserData<number, dim> >, SmartPtr<CplUserData<number, dim> >)>(&T::set_conductances), "", "conductance functions for K+, Na+ and leakage (in S/m^2)")
			#ifdef UG_FOR_LUA
				.add_method("set_conductances", static_cast<void (T::*)(const char*, const char* , const char*)>(&T::set_conductances), "", "conductance function names for K+, Na+ and leakage (in S/m^2)")
			#endif

				.add_method("set_rev_pot", static_cast<void (T::*)(number, number, number)>(&T::set_rev_pot), "", "reversal potentials (in V) for K+, Na+ and leakage")
				.add_method("set_rev_pot", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, SmartPtr<CplUserData<number, dim> >, SmartPtr<CplUserData<number, dim> >)>(&T::set_rev_pot), "", "reversal potential (in V) functions for K+, Na+ and leakage")
			#ifdef UG_FOR_LUA
				.add_method("set_rev_pot", static_cast<void (T::*)(const char*, const char* , const char*)>(&T::set_rev_pot), "", "reversal potential (in V) function names for K+, Na+ and leakage")
			#endif

				.add_method("set_diameter", &T::set_diameter, "", "constant diameter (in m)")
				.add_method("set_injection", &T::set_injection, "", "injection current (in A/m^2)")
				.add_method("set_temperature", &T::set_temperature, "", "temperature (in K)")
				.add_method("enable_temperature_dependency", &T::enable_temperature_dependency,
					"", "whether gating change rates are to be temperature-dependent");
			reg.add_class_to_group(name, "ElemDiscHH_Base", tag);
		}

		// fully implicit cable discretization with constant reversal potential for K and Na
		{
			typedef ImplicitActiveCableDisc<TDomain> T;
			typedef ImplicitActiveCableDiscBase<TDomain> TBase;
			string name = string("ImplicitActiveCableDisc").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*) (const char*, const char*)>("function(s) # subset(s)")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ImplicitActiveCableDisc", tag);
		}

		// Kabel Diff FV1 with dynamically calculated K+ and Na+ Nernst potentials
		{
			typedef ImplicitActiveCableDiscNernst<TDomain> T;
			typedef ImplicitActiveCableDiscBase<TDomain> TBase;
			string name = string("ImplicitActiveCableDiscNernst").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("function(s) # subset(s)")
				.add_method("set_diffusion_constants", &T::set_diffusion_constants)
				.add_method("set_outside_concs", &T::set_outside_concs)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ImplicitActiveCableDiscNernst", tag);
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
				.add_method("enable_temperature_dependency", &T::enable_temperature_dependency)
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
							"sets leak conductance for leakage")
				.add_method("set_cond", static_cast<void (T::*)(number, const char*)>(&T::set_cond),
							"", "leak conductance (S/m^2) | default | value=1.0 # subset(s) as C-type string",
							"sets leak conductance for leakage")
				.add_method("set_cond", static_cast<void (T::*)(number, const std::vector<std::string>&)>(&T::set_cond),
							"", "leak conductance (S/m^2) | default | value=1.0 # subset(s) as vector of strings",
							"sets leak conductance for leakage")
#ifdef UG_FOR_LUA
				.add_method("set_cond", static_cast<void (T::*) (SmartPtr<LuaUserData<number, dim> >)>(&T::set_cond),
					"", "leak conductance function (S/m^2)", "sets coordinate-dependent leak conductance value")
				.add_method("set_cond", static_cast<void (T::*) (const char*)>(&T::set_cond),
					"", "Lua name of leak conductance function (S/m^2)",
					"sets coordinate-dependent leak conductance value")
#endif
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
				.add_method("set_perm", static_cast<void (T::*)(number)>(&T::set_perm),
					"", "permeability", "sets a constant permeability")
				.add_method("set_perm", static_cast<void (T::*)(number, number, number, number, int)>(&T::set_perm),
					"", "flux at rest (mol/(m^2*s))#inner concentration at rest (mM)# outer concentration at rest (mM)#potential at rest (V)",
					"tunes the channel to equilibrate fluxes at resting conditions")
				.add_method("set_valency", &T::set_valency, "", "ion valency", "")
				.add_method("set_ohmic", &T::set_ohmic, "", "ohmic?", "whether an Ohmic model is to be used")
				.add_method("set_cond", &T::set_cond, "", "conductance", "set the conductance")
				.add_method("set_rev_pot", &T::set_rev_pot, "", "reversal potential", "set the reversal potential")
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


		// K-channel (A-type) Golding 2001
		{
			typedef KA_Golding01<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("KA_Golding01").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductance", static_cast<void (T::*) (number)>(&T::set_gkbar),
					"", "K conductance (S/m^2)", "")
#ifdef UG_FOR_LUA
				.add_method("set_conductance", static_cast<void (T::*) (SmartPtr<LuaUserData<number, dim> >)>(&T::set_gkbar),
					"", "K conductance function (S/m^2)", "sets coordinate-dependent K conductance value")
				.add_method("set_conductance", static_cast<void (T::*) (const char*)>(&T::set_gkbar),
					"", "Lua name of K conductance function (S/m^2)",
					"sets coordinate-dependent K conductance value")
#endif
				.add_method("set_vhalfn", &T::set_vhalfn)
				.add_method("set_a0n", &T::set_a0n)
				.add_method("set_zetan", &T::set_zetan)
				.add_method("set_gmn", &T::set_gmn)
				.add_method("set_vhalfn_dist", &T::set_vhalfn_dist)
				.add_method("set_a0n_dist", &T::set_a0n_dist)
				.add_method("set_zetan_dist", &T::set_zetan_dist)
				.add_method("set_gmn_dist", &T::set_gmn_dist)
#ifdef UG_FOR_LUA
				// investigate: very strangely, it is not possible to register the following two methods in reverse order...:
				// this gives: Registry ERROR: Unregistered Class used in Method:
				// 'void KA_Golding011d:set_proximality_fct(SmartPtr<> set proximality function)': for Parameter 1
				.add_method("set_proximality_fct", static_cast<void (T::*) (const char*)>(&T::set_proximality_fct),
					"", "proximality function name",
					"The proximality function returns 1 if a position is to be considered proximal "
					"and 0 if it is to be considered distal.")
				.add_method("set_proximality_fct", static_cast<void (T::*) (SmartPtr<LuaUserData<number, dim> >)>(&T::set_proximality_fct),
					"", "proximality function",
					"The proximality function returns 1 if a position is to be considered proximal "
					"and 0 if it is to be considered distal.")
#endif
				.add_method("set_logNGate", &T::set_logNGate)
				.add_method("set_logLGate", &T::set_logLGate)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "KA_Golding01", tag);
		}

		// K-channel (delayed rectifier) Golding 2001
		{
			typedef KDR_Golding01<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("KDR_Golding01").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductance", static_cast<void (T::*) (number)>(&T::set_gkbar), "",
					"K conductance (S/m^2) | default | value=30.0#"
					"subset(s) as vector of string", "sets K conductance value")
#ifdef UG_FOR_LUA
				.add_method("set_conductance", static_cast<void (T::*) (SmartPtr<LuaUserData<number, dim> >)>(&T::set_gkbar),
					"", "K conductance function (S/m^2)", "sets coordinate-dependent K conductance value")
				.add_method("set_conductance", static_cast<void (T::*) (const char*)>(&T::set_gkbar),
					"", "Lua name of K conductance function (S/m^2)",
					"sets coordinate-dependent K conductance value")
#endif
				.add_method("set_logNGate", &T::set_logNGate)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "KDR_Golding01", tag);
		}

		// Na-channel (axon) Golding 2001
		{
			typedef Nax_Golding01<TDomain> T;
			typedef ICableMembraneTransport<TDomain> TBase;
			string name = string("Nax_Golding01").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductance", static_cast<void (T::*) (number)>(&T::set_gbar), "",
					"Na conductance (S/m^2) | default | value=100.0#"
					"subset(s) as vector of string", "sets Na conductance value")
#ifdef UG_FOR_LUA
				.add_method("set_conductance", static_cast<void (T::*) (SmartPtr<LuaUserData<number, dim> >)>(&T::set_gbar),
					"", "Na conductance function (S/m^2)", "sets coordinate-dependent Na conductance value")
				.add_method("set_conductance", static_cast<void (T::*) (const char*)>(&T::set_gbar),
					"", "Lua name of Na conductance function (S/m^2)",
					"sets coordinate-dependent Na conductance value")
#endif
				.add_method("set_logMGate", &T::set_logMGate)
				.add_method("set_logHGate", &T::set_logHGate)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Nax_Golding01", tag);
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
#ifdef UG_FOR_LUA
				.add_method("set_spec_cap", static_cast<void (T::*)(SmartPtr<LuaUserData<number, TDomain::dim> >)>(&T::set_spec_cap),
					"", "specific capacitance function (F/m^2)",
					"sets a new specific capacitance on subsets given as C-style string")
				.add_method("set_spec_cap", static_cast<void (T::*)(const char*)>(&T::set_spec_cap),
					"", "specific capacitance function name (F/m^2)",
					"sets a new specific capacitance on subsets given as vector of string")
#endif

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
					"", "temperature (degrees C) | default | value=37", "sets new temperature")

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
#ifdef UG_FOR_LUA
				.add_method("set_influx_function",
					static_cast<void (T::*)(const char*, const char*)> (&T::set_influx_function),
					"", "lua function name # subset name", "Set an influx density function on a subset.")
#ifndef UG_FOR_VRL
				.add_method("set_influx_function",
					static_cast<void (T::*)(SmartPtr<LuaUserData<number, dim> >, const std::string&)> (&T::set_influx_function),
					"", "LuaUserFunction # subset name", "Set an influx density function on a subset.")
#endif
#endif
				.add_method("set_influx_subset", &T::set_influx_subset)
#ifdef UG_CPU_1
				.add_method("estimate_cfl_cond", &T::template estimate_cfl_cond<CPUAlgebra::vector_type>)
#endif
				.add_method("set_synapse_handler", &T::set_synapse_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "CableEquation", tag);
		}

		// CableEquationWithOuterPot
		{
			typedef CableEquationWithOuterPot<TDomain> T;
			typedef CableEquation<TDomain> TBase;
			string name = string("CableEquationWithOuterPot").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>
					("Fcts (inner potential, outer potential) # Subset(s)")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "CableEquationWithOuterPot", tag);
		}

#ifdef CONVERTED_TRANSPORT_ENABLED
		#include "membrane_transport/nmodl_converter/converted/channels.cpp"
#endif

		// Cuthill-McKee ordering
		{
			reg.add_function("order_cuthillmckee", &order_cuthillmckee<TDomain>, grp.c_str(),
							 "", "approxSpace", "vertex ordering for solver optimization");
		}

		// domain scaling
		{
			reg.add_function("scale_domain", &scale_domain<TDomain>, grp.c_str(),
							 "", "domain, scaling factor", "scales the domain (also diameter)");
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
			//reg.add_function("test_vertices", &test_vertices<TDomain>, grp.c_str(),
			//				 "", "domain", "Tests the distributed vertices for correctness of the virtual table pointer.");
		}

		{
			typedef SynapseHandler<TDomain> TSH;
			string name = string("SynapseHandler").append(suffix);
			reg.add_class_<TSH>(name, grp)
				.template add_constructor<void (*)()>()
                .add_method("set_ce_object", &TSH::set_ce_object)

                .add_method("set_activation_timing_alpha",
                    static_cast<void (TSH::*)(number, number, number, number)>(&TSH::set_activation_timing_alpha),
                    "", "mean onset#mean tau#onset deviation#tau deviation", "")
                .add_method("set_activation_timing_alpha",
                    static_cast<void (TSH::*)(number, number, number, number, number)>(&TSH::set_activation_timing_alpha),
                    "", "mean onset#mean tau#onset deviation#tau deviation#peak conductance", "")
                .add_method("set_activation_timing_alpha",
                    static_cast<void (TSH::*)(number, number, number, number, number, bool)>(&TSH::set_activation_timing_alpha),
                    "", "mean onset#mean tau#onset deviation#tau deviation#peak conductance#constant seed?", "")
                .add_method("set_activation_timing_alpha",
                    static_cast<void (TSH::*)(number, number, number, number, number, number, bool)>(&TSH::set_activation_timing_alpha),
                    "", "mean onset#mean tau#mean peak conductance#onset deviation#tau deviation#peak conductance deviation#constant seed?", "")

                .add_method("set_activation_timing_biexp",
                    static_cast<void (TSH::*)(number, number, number, number, number, number)>(&TSH::set_activation_timing_biexp),
                    "", "mean onset#mean tau1#mean tau2#onset deviation#tau1 deviation#tau2 deviation", "")
                .add_method("set_activation_timing_biexp",
                    static_cast<void (TSH::*)(number, number, number, number, number, number, number)>(&TSH::set_activation_timing_biexp),
                    "", "mean onset#mean tau1#mean tau2#onset deviation#tau1 deviation#tau2 deviation#peak conductance", "")
                .add_method("set_activation_timing_biexp",
                    static_cast<void (TSH::*)(number, number, number, number, number, number, number, bool)>(&TSH::set_activation_timing_biexp),
                    "", "mean onset#mean tau1#mean tau2#onset deviation#tau1 deviation#tau2 deviation#peak conductance#constant seed?", "")
                .add_method("set_activation_timing_biexp",
                    static_cast<void (TSH::*)(number, number, number, number, number, number, number, number, bool)>(&TSH::set_activation_timing_biexp),
                    "", "mean onset#mean tau1#mean tau2#mean peak conductance#onset deviation#tau1 deviation#tau2 deviation#peak conductance deviation#constant seed?", "")

                .add_method("add_activation_timing_alpha_ball", &TSH::add_activation_timing_alpha_ball,
                    "", "timing values as six-component vector (mean onset, dev onset, mean tau, dev tau, mean peak conductance, dev peak conductance)"
                        "#ball region as four component vector (center coordinates, radius)",
                    "Add a ball-shaped region with specific activation pattern for alpha post-synapses.")
                 .add_method("add_activation_timing_exp2_ball", &TSH::add_activation_timing_exp2_ball,
                                        "", "timing values as six-component vector (mean onset, mean tau1, mean tau2, mean peak conductance, onset dev, tau1 dev, tau2 dev, dev peak conductance)"
                                            "#ball region as four component vector (center coordinates, radius)",
                                        "Add a ball-shaped region with specific activation pattern for exp2 post-synapses.")

				.add_method("show_status", &TSH::show_status)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "SynapseHandler", tag);
		}


#ifdef CONVERTED_TRANSPORT_BAP_ENABLED
		{
			typedef nax_g01_converted_standard_UG<TDomain> TNAX;
			string name = string("nax_g01").append(suffix);
			reg.add_class_<TNAX>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>()
				.add_method("init", &TNAX::init)
				.add_method("init_attachments", &TNAX::init_attachments)
				.add_method("update_gating", &TNAX::update_gating)
				.add_method("current", &TNAX::current)
				.add_method("ce_obj_available", &TNAX::ce_obj_available)
				.add_method("state_values", &TNAX::state_values)
				.add_method("getgbar", &TNAX::getgbar)
				.add_method("gettha", &TNAX::gettha)
				.add_method("getqa", &TNAX::getqa)
				.add_method("getRa", &TNAX::getRa)
				.add_method("getRb", &TNAX::getRb)
				.add_method("getthi1", &TNAX::getthi1)
				.add_method("getthi2", &TNAX::getthi2)
				.add_method("getqd", &TNAX::getqd)
				.add_method("getqg", &TNAX::getqg)
				.add_method("getmmin", &TNAX::getmmin)
				.add_method("gethmin", &TNAX::gethmin)
				.add_method("getq10", &TNAX::getq10)
				.add_method("getRg", &TNAX::getRg)
				.add_method("getRd", &TNAX::getRd)
				.add_method("getthinf", &TNAX::getthinf)
				.add_method("getqinf", &TNAX::getqinf)
				.add_method("getena", &TNAX::getena)
				.add_method("getmscale", &TNAX::getmscale)
				.add_method("gethscale", &TNAX::gethscale)
				.add_method("setgbar", &TNAX::setgbar)
				.add_method("settha", &TNAX::settha)
				.add_method("setqa", &TNAX::setqa)
				.add_method("setRa", &TNAX::setRa)
				.add_method("setRb", &TNAX::setRb)
				.add_method("setthi1", &TNAX::setthi1)
				.add_method("setthi2", &TNAX::setthi2)
				.add_method("setqd", &TNAX::setqd)
				.add_method("setqg", &TNAX::setqg)
				.add_method("setmmin", &TNAX::setmmin)
				.add_method("sethmin", &TNAX::sethmin)
				.add_method("setq10", &TNAX::setq10)
				.add_method("setRg", &TNAX::setRg)
				.add_method("setRd", &TNAX::setRd)
				.add_method("setthinf", &TNAX::setthinf)
				.add_method("setqinf", &TNAX::setqinf)
				.add_method("setena", &TNAX::setena)
				.add_method("setmscale", &TNAX::setmscale)
				.add_method("sethscale", &TNAX::sethscale)
				.add_method("set_log_mGate", &TNAX::set_log_mGate)
				.add_method("set_log_hGate", &TNAX::set_log_hGate)

				;
		}
		{
			typedef h_g05_converted_standard_UG<TDomain> TH;
			string name = string("h_g05").append(suffix);
			reg.add_class_<TH>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>()
				.add_method("alpha", &TH::alpha)
				.add_method("beta", &TH::beta)
				.add_method("init", &TH::init)
				.add_method("update_gating", &TH::update_gating)
				.add_method("current", &TH::current)
				.add_method("ce_obj_available", &TH::ce_obj_available)
				.add_method("state_values", &TH::state_values)
				.add_method("getgbar", &TH::getgbar)
				.add_method("geterevh", &TH::geterevh)
				.add_method("getvhalf", &TH::getvhalf)
				.add_method("geta0", &TH::geta0)
				.add_method("getzeta", &TH::getzeta)
				.add_method("getab", &TH::getab)
				.add_method("getqten", &TH::getqten)
				.add_method("gettemp", &TH::gettemp)
				.add_method("getgas", &TH::getgas)
				.add_method("getfarad", &TH::getfarad)
				.add_method("setgbar", &TH::setgbar)
				.add_method("seterevh", &TH::seterevh)
				.add_method("setvhalf", &TH::setvhalf)
				.add_method("seta0", &TH::seta0)
				.add_method("setzeta", &TH::setzeta)
				.add_method("setab", &TH::setab)
				.add_method("setqten", &TH::setqten)
				.add_method("settemp", &TH::settemp)
				.add_method("setgas", &TH::setgas)
				.add_method("setfarad", &TH::setfarad)
				.add_method("set_log_SGate", &TH::set_log_SGate)
				.add_method("set_log_hhGate", &TH::set_log_hhGate)
					;

		}
		{
			typedef kadist_g01_converted_standard_UG<TDomain> TH;
			string name = string("kadist_g01").append(suffix);
			reg.add_class_<TH>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>()
				.add_method("init", &TH::init)
				.add_method("update_gating", &TH::update_gating)
				.add_method("current", &TH::current)
				.add_method("ce_obj_available", &TH::ce_obj_available)
				.add_method("state_values", &TH::state_values)
				.add_method("alpn", &TH::alpn)
				.add_method("betn", &TH::betn)
				.add_method("alpl", &TH::alpl)
				.add_method("betl", &TH::betl)
				.add_method("getek", &TH::getek)
				.add_method("getgkabar", &TH::getgkabar)
				.add_method("getvhalfn", &TH::getvhalfn)
				.add_method("getvhalfl", &TH::getvhalfl)
				.add_method("geta0l", &TH::geta0l)
				.add_method("geta0n", &TH::geta0n)
				.add_method("getzetan", &TH::getzetan)
				.add_method("getzetal", &TH::getzetal)
				.add_method("getgmn", &TH::getgmn)
				.add_method("getgml", &TH::getgml)
				.add_method("getlmin", &TH::getlmin)
				.add_method("getnmin", &TH::getnmin)
				.add_method("getpw", &TH::getpw)
				.add_method("gettq", &TH::gettq)
				.add_method("getqq", &TH::getqq)
				.add_method("getq10", &TH::getq10)
				.add_method("getqtl", &TH::getqtl)
				.add_method("getnscale", &TH::getnscale)
				.add_method("getlscale", &TH::getlscale)
				.add_method("setek", &TH::setek)
				.add_method("setgkabar", &TH::setgkabar)
				.add_method("setvhalfn", &TH::setvhalfn)
				.add_method("setvhalfl", &TH::setvhalfl)
				.add_method("seta0l", &TH::seta0l)
				.add_method("seta0n", &TH::seta0n)
				.add_method("setzetan", &TH::setzetan)
				.add_method("setzetal", &TH::setzetal)
				.add_method("setgmn", &TH::setgmn)
				.add_method("setgml", &TH::setgml)
				.add_method("setlmin", &TH::setlmin)
				.add_method("setnmin", &TH::setnmin)
				.add_method("setpw", &TH::setpw)
				.add_method("settq", &TH::settq)
				.add_method("setqq", &TH::setqq)
				.add_method("setq10", &TH::setq10)
				.add_method("setqtl", &TH::setqtl)
				.add_method("setnscale", &TH::setnscale)
				.add_method("setlscale", &TH::setlscale)
				.add_method("set_log_SGate", &TH::set_log_SGate)
				.add_method("set_log_nGate", &TH::set_log_nGate)
				.add_method("set_log_lGate", &TH::set_log_lGate)
				;
		}
		{
			typedef kaprox_g01_converted_standard_UG<TDomain> TH;
			string name = string("kaprox_g01").append(suffix);
			reg.add_class_<TH>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>()
				.add_method("init", &TH::init)
				.add_method("update_gating", &TH::update_gating)
				.add_method("current", &TH::current)
				.add_method("ce_obj_available", &TH::ce_obj_available)
				.add_method("state_values", &TH::state_values)
				.add_method("alpn", &TH::alpn)
				.add_method("betn", &TH::betn)
				.add_method("alpl", &TH::alpl)
				.add_method("betl", &TH::betl)
				.add_method("getek", &TH::getek)
				.add_method("getgkabar", &TH::getgkabar)
				.add_method("getvhalfn", &TH::getvhalfn)
				.add_method("getvhalfl", &TH::getvhalfl)
				.add_method("geta0l", &TH::geta0l)
				.add_method("geta0n", &TH::geta0n)
				.add_method("getzetan", &TH::getzetan)
				.add_method("getzetal", &TH::getzetal)
				.add_method("getgmn", &TH::getgmn)
				.add_method("getgml", &TH::getgml)
				.add_method("getlmin", &TH::getlmin)
				.add_method("getnmin", &TH::getnmin)
				.add_method("getpw", &TH::getpw)
				.add_method("gettq", &TH::gettq)
				.add_method("getqq", &TH::getqq)
				.add_method("getq10", &TH::getq10)
				.add_method("getqtl", &TH::getqtl)
				.add_method("getnscale", &TH::getnscale)
				.add_method("getlscale", &TH::getlscale)
				.add_method("setek", &TH::setek)
				.add_method("setgkabar", &TH::setgkabar)
				.add_method("setvhalfn", &TH::setvhalfn)
				.add_method("setvhalfl", &TH::setvhalfl)
				.add_method("seta0l", &TH::seta0l)
				.add_method("seta0n", &TH::seta0n)
				.add_method("setzetan", &TH::setzetan)
				.add_method("setzetal", &TH::setzetal)
				.add_method("setgmn", &TH::setgmn)
				.add_method("setgml", &TH::setgml)
				.add_method("setlmin", &TH::setlmin)
				.add_method("setnmin", &TH::setnmin)
				.add_method("setpw", &TH::setpw)
				.add_method("settq", &TH::settq)
				.add_method("setqq", &TH::setqq)
				.add_method("setq10", &TH::setq10)
				.add_method("setqtl", &TH::setqtl)
				.add_method("setnscale", &TH::setnscale)
				.add_method("setlscale", &TH::setlscale)
				.add_method("set_log_SGate", &TH::set_log_SGate)
				.add_method("set_log_nGate", &TH::set_log_nGate)
				.add_method("set_log_lGate", &TH::set_log_lGate)
				;
		}
		{
			typedef kdrca1_g01_converted_standard_UG<TDomain> TH;
			string name = string("kdrca1_g01").append(suffix);
			reg.add_class_<TH>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>()
				.add_method("init", &TH::init)
				.add_method("update_gating", &TH::update_gating)
				.add_method("current", &TH::current)
				.add_method("ce_obj_available", &TH::ce_obj_available)
				.add_method("state_values", &TH::state_values)
				.add_method("alpn", &TH::alpn)
				.add_method("betn", &TH::betn)
				.add_method("getek", &TH::getek)
				.add_method("getgkdrbar", &TH::getgkdrbar)
				.add_method("getikmax", &TH::getikmax)
				.add_method("getvhalfn", &TH::getvhalfn)
				.add_method("geta0n", &TH::geta0n)
				.add_method("getzetan", &TH::getzetan)
				.add_method("getgmn", &TH::getgmn)
				.add_method("getnmax", &TH::getnmax)
				.add_method("getq10", &TH::getq10)
				.add_method("getnscale", &TH::getnscale)
				.add_method("setek", &TH::setek)
				.add_method("setgkdrbar", &TH::setgkdrbar)
				.add_method("setikmax", &TH::setikmax)
				.add_method("setvhalfn", &TH::setvhalfn)
				.add_method("seta0n", &TH::seta0n)
				.add_method("setzetan", &TH::setzetan)
				.add_method("setgmn", &TH::setgmn)
				.add_method("setnmax", &TH::setnmax)
				.add_method("setq10", &TH::setq10)
				.add_method("setnscale", &TH::setnscale)
				.add_method("set_log_SGate", &TH::set_log_SGate)
				.add_method("set_log_nGate", &TH::set_log_nGate)
					;
		}
		{
			//todo: spines?
		}
#endif
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
        {
            typedef synapse_handler::IBaseSynapse TBS;
            reg.add_class_<TBS>("IBaseSynapse", grp)
                    ;
        }
        {
            typedef synapse_handler::IPreSynapse T;
            typedef synapse_handler::IBaseSynapse TBase;
            reg.add_class_<T, TBase>("IPreSynapse", grp);
        }
        {
            typedef synapse_handler::IPostSynapse T;
            typedef synapse_handler::IBaseSynapse TBase;
            reg.add_class_<T, TBase>("IPostSynapse", grp);
        }

		{
			string name = "ISynapseContainer";

			typedef ISynapseContainer TSC;
			reg.add_class_<TSC>(name, grp);;
		}

        {
            string name = "AlphaSynapses";

            typedef AlphaSynapses TSC;
            reg.add_class_<TSC>(name, grp)
                .add_constructor<void (*)(const size_t, const size_t)>("infile#outfile", "Initializes a SynapseDistributor", grp, "")
                //.add_constructor<void (*)(const size_t, number, number, number, number, number, number, number)>("dom#outfile", "Initializes a SynapseDistributor", grp, "")
                .add_method("set_mean_gMax",&TSC::set_mean_gMax)
                .add_method("set_dev_gMax",&TSC::set_dev_gMax)
                .add_method("set_mean_onset",&TSC::set_mean_onset)
                .add_method("set_dev_onset",&TSC::set_dev_onset)
                .add_method("set_mean_tau",&TSC::set_mean_tau)
                .add_method("set_dev_tau",&TSC::set_dev_tau)
                .add_method("set_mean_e",&TSC::set_mean_rev)
                .add_method("set_dev_e",&TSC::set_dev_rev)
                .add_method("size",&TSC::size)
                .add_method("get_synapses",&TSC::get_synapses);
        }

		{
			string name = "Exp2Synapses";

			typedef Exp2Synapses TSC;
			reg.add_class_<TSC>(name, grp)
				.add_constructor<void (*)(const size_t, const size_t)>("infile#outfile", "Initializes a SynapseDistributor", grp, "")
				//.add_constructor<void (*)(const size_t, number, number, number, number, number, number, number)>("dom#outfile", "Initializes a SynapseDistributor", grp, "")
				.add_method("set_mean_onset",&TSC::set_mean_onset)
				.add_method("set_dev_onset",&TSC::set_dev_onset)
				.add_method("set_mean_tau1",&TSC::set_mean_tau1)
				.add_method("set_dev_tau1",&TSC::set_dev_tau1)
				.add_method("set_mean_tau2",&TSC::set_mean_tau2)
				.add_method("set_dev_tau2",&TSC::set_dev_tau2)
				.add_method("set_mean_e",&TSC::set_mean_rev)
				.add_method("set_dev_e",&TSC::set_dev_rev)
				.add_method("set_mean_w",&TSC::set_mean_gMax)
				.add_method("set_dev_w",&TSC::set_dev_gMax)
				.add_method("size",&TSC::size)
				.add_method("get_synapses",&TSC::get_synapses);
		}

        {
            string name = "AlphaSynapsePair";

            typedef AlphaSynapsePair T;
            reg.add_class_<T>(name, grp)
                .add_constructor()
                .add_method("set_id",&T::set_id)
                .add_method("set_onset",&T::set_onset)
                .add_method("set_tau",&T::set_tau)
                .add_method("set_gMax",&T::set_gMax)
                .add_method("set_reversal_potential",&T::set_reversal_potential)
                .add_method("pre_synapse",&T::pre_synapse)
                .add_method("post_synapse",&T::post_synapse)
                .set_construct_as_smart_pointer(true);
        }

        {
            string name = "Exp2SynapsePair";

            typedef Exp2SynapsePair T;
            reg.add_class_<T>(name, grp)
                .add_constructor()
                .add_method("set_id",&T::set_id)
                .add_method("set_threshold",&T::set_threshold)
                .add_method("set_gMax",&T::set_gMax)
                .add_method("set_reversal_potential",&T::set_reversal_potential)
                .add_method("set_taus",&T::set_taus)
                .add_method("pre_synapse",&T::pre_synapse)
                .add_method("post_synapse",&T::post_synapse)
                .set_construct_as_smart_pointer(true);
        }


		{
			// ////////////////////////////////////////////
			// //////// SynapseDistributor ////////////////
			// ////////////////////////////////////////////

			// TODO: might better be registered (and implemented) in a domain-dependent manner
			string name = "SynapseDistributor";
			string constr_params("infile#outfile");

			typedef SynapseDistributor TSD;

			reg.add_class_<TSD>(name, grp)
				.add_constructor<void (*)(string)>("infile")
				.add_constructor<void (*)(SmartPtr<ug::Domain3d>)>("dom")

				.add_method("clear", static_cast<void (TSD::*)()>(&TSD::clear), "", "",
				    "Removes all synapses from grid.")
				.add_method("clear", static_cast<void (TSD::*)(int)>(&TSD::clear), "", "subsetIndex",
				    "Removes all Synapses from subset.")

				.add_method("place_synapse_at_coords", &TSD::place_synapse_at_coords, "",
				    "coordinates (as three-component vector)#pre-synapse#post-synapse",
				    "Places one pair of pre- and post-synapse on the edge nearest to the given coordinates.")

				.add_method("place_synapses_uniform",
                    static_cast<void (TSD::*)(size_t, const string&)>(&TSD::place_synapses_uniform),
                    "", "number of synapses#synapse type", "Distributes post-synapses uniformly on grid.")
                .add_method("place_synapses_uniform",
                    static_cast<void (TSD::*)(int, size_t, const string&)>(&TSD::place_synapses_uniform),
                    "", "subset index#number of synapses#synapse type", "Distributes post-synapses uniformly on subset.")
                .add_method("place_synapses_uniform",
                    static_cast<void (TSD::*)(const char*, size_t, const string&)>(&TSD::place_synapses_uniform),
                    "", "subset name#number of synapses#synapse type", "Distributes post-synapses uniformly on subset.")
                .add_method("place_synapses_uniform_density",
                    static_cast<void (TSD::*)(int, number, const string&)>(&TSD::place_synapses_uniform),
                    "", "subset index#density (m^-1)#synapse type", "Uniformly distributes post-synapses on subset with given density.")
                .add_method("place_synapses_uniform",
                    static_cast<void (TSD::*)(const char*, number, const string&)>(&TSD::place_synapses_uniform),
                    "", "subset name#density (m^-1)#synapse type", "Uniformly distributes post-synapses on subset with given density.")
				.add_method("place_synapses_uniform_density",
                    static_cast<void (TSD::*)(number, number, number, number, number, const string&)>(&TSD::place_synapses_uniform),
                    "", "subset name#density (m^-1)#synapse type", "Uniformly distributes post-synapses in ball region.")
                .add_method("place_synapses_uniform",
                    static_cast<void (TSD::*)(size_t, number, number, number, number, const string&)>(&TSD::place_synapses_uniform),
                    "", "number of synapses#x#y#z#radius#type of synapse")
				.add_method("place_synapses", &TSD::place_synapses,
				    "", "", "Distributes post-synapses according to given distribution on the subsets.")
				.add_method("degenerate_uniform", static_cast<void (TSD::*)(number)>(&TSD::degenerate_uniform),
				    "", "percentage", "Removes a percentage of synapses from the grid.")
				.add_method("degenerate_uniform", static_cast<void (TSD::*)(number, int)>(&TSD::degenerate_uniform),
				    "", "percentage#subset index", "Removes a percentage of synapses from the given subset.")
				.add_method("degenerate_uniform", static_cast<void (TSD::*)(number, const char*)>(&TSD::degenerate_uniform),
				    "", "percentage#subset name", "Removes a percentage of synapses from the given subset.",grp)

				.add_method("num_synapses", static_cast<size_t (TSD::*)() const>(&TSD::num_synapses),
				    "", "", "Returns global number of synapses.")
				.add_method("num_synapses", static_cast<size_t (TSD::*)(int) const>(&TSD::num_synapses),
				    "", "subset index", "Returns number of synapses in specified subset.")
				.add_method("num_synapses", static_cast<size_t (TSD::*)(const char*) const>(&TSD::num_synapses),
				    "", "subset name", "Returns number of synapses in specified subset.")

				.add_method("print_status", &TSD::print_status, "", "", "Prints synapse status of grid.")
				.add_method("export_grid", &TSD::export_grid, "", "file name", "Saves grid with synapses to file.")

				.set_construct_as_smart_pointer(true);
		}

		// ugx -> swc conversion
		reg.add_function("save_neuron_to_swc", &save_neuron_to_swc, grp.c_str(),
		                 "", "ugx file name # neuron index # swc file name # scale",
		                 "Save a single neuron from a ugx neuronal network to swc.");

		// innermost_neuron_id_in_subset
		reg.add_function("innermost_neuron_id_in_subset", &innermost_neuron_id_in_subset, grp.c_str(),
						 "", "subset name # subset handler", "");

		// neurite subset length
		reg.add_function("subset_length", static_cast<number (*) (int, ConstSmartPtr<MGSubsetHandler>)>(&subset_length),
            grp.c_str(), "length of neurite subset",
            "subset index#subset handler", "Get the total length of all neurite segments in subset.");
		reg.add_function("subset_length", static_cast<number (*) (const char*, ConstSmartPtr<MGSubsetHandler>)>(&subset_length),
            grp.c_str(), "length of neurite subset",
            "subset index#subset handler", "Get the total length of all neurite segments in subset.");

#ifdef UG_PARALLEL
#ifdef UG_DIM_3
		{
			typedef DomainPartitioner<Domain3d, SimpleCableNeuronPartitioner> T;
			string name = "SimpleCableNeuronPartitioner";
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(Domain3d&)>("domain")
				.set_construct_as_smart_pointer(true);
		}
#endif
#endif

#ifdef NC_WITH_PARMETIS
		{
			typedef CableNeuronUnificator T;
			std::string name = std::string("CableNeuronUnificator");
			typedef parmetis::IUnificator<Edge> TBase;
			reg.add_class_<T, TBase>(name, grp)
				.add_constructor<void (*)()>()
				.set_construct_as_smart_pointer(true);
		}
#endif

		// NeuronalTopologyImporter
		{
			typedef neuronal_topology_importer::NeuronalTopologyImporter T;
			std::string name = std::string("NeuronalTopologyImporter");
			reg.add_class_<T>(name, grp)
				.add_constructor<void (*) ()>("", "", grp)
				.add_method("set_scaling_factors", &T::set_scaling_factors, "", "length # time # potential # conductance", "")
				.add_method("set_length_scaling", &T::set_length_scaling, "", "length scaling factor", "")
				.add_method("set_time_scaling", &T::set_time_scaling, "", "time scaling factor", "")
				.add_method("set_potential_scaling", &T::set_potential_scaling, "", "potential scaling factor", "")
				.add_method("set_conductance_scaling", &T::set_conductance_scaling, "", "conductance scaling factor", "")
				.add_method("add_joining_criterion", &T::add_joining_criterion, "", "criterion", "")
				.add_method("rm_joining_criterion", &T::rm_joining_criterion, "", "criterion", "")
				.add_method("import_geometry_and_generate_grid",
					(bool (T::*)(const std::string&)) (&T::import_and_generate_grid),
					"success", "file name")
				.add_method("import_geometry_and_generate_grid",
					(bool (T::*)(const std::string&, const std::string&)) (&T::import_and_generate_grid),
					"success", "file base name # format (hoc, ngx, swc, txt)")
				.set_construct_as_smart_pointer(true);
		}
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
	typedef Attachment<std::vector<IBaseSynapse*> >AVSynapse;
	typedef Attachment<uint> ANeuronID;

	GlobalAttachments::declare_attachment<ADiameter>("diameter", true);
	GlobalAttachments::declare_attachment<AVSynapse>("synapses", false);
	GlobalAttachments::declare_attachment<ANeuronID>("neuronID", false);


	typedef cable_neuron::Functionality Functionality;

	// we only need to compile for 3D, as cable models are always in a 3D space
	typedef boost::mpl::list<
	#ifdef UG_DIM_3
		Domain3d
	#endif
	> OurCompileDomainList;

	try
	{
		RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality, OurCompileDomainList>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality, OurCompileDomainList>(*reg,grp);
	}
    UG_REGISTRY_CATCH_THROW(grp);


    // Register synapse types.
    // This has to be done only AFTER the rest has been registered.
    // Synapse iterators and begin() and end() methods of SynapseHandler
    // are automatically registered as well.
	try
	{
	    {
            typedef cable_neuron::synapse_handler::OnsetPreSynapse T;
            SynapseDealer::instance()->register_pre_synapse_type<T>(reg, grp)
                .add_method("onset", &T::onset)
                .add_method("duration", &T::duration)
                .add_method("set_onset", &T::set_onset)
                .add_method("set_duration", &T::set_duration)
                .add_method("name", &T::name);
	    }
	    {
	        typedef cable_neuron::synapse_handler::ThresholdPreSynapse T;
	        SynapseDealer::instance()->register_pre_synapse_type<T>(reg, grp)
                .add_method("set_duration", &T::set_duration)
                .add_method("set_threshold", &T::set_threshold)
                .add_method("name", &T::name);
	    }
	    {
            typedef cable_neuron::synapse_handler::AlphaPostSynapse T;
	        SynapseDealer::instance()->register_post_synapse_type<T>(reg, grp)
                .add_method("tau", &T::tau)
                .add_method("gMax", &T::gMax)
                .add_method("rev", &T::rev)
                .add_method("set_tau", &T::set_tau)
                .add_method("set_gMax", &T::set_gMax)
                .add_method("set_rev", &T::set_rev)
                .add_method("name", &T::name);
	    }
	    {
            typedef cable_neuron::synapse_handler::Exp2PostSynapse T;
	        SynapseDealer::instance()->register_post_synapse_type<T>(reg, grp)
                .add_method("rev", &T::rev)
                .add_method("gMax", &T::gMax)
                .add_method("tau1", &T::tau1)
                .add_method("tau2", &T::tau2)
                .add_method("set_rev", &T::set_rev)
                .add_method("set_gMax", &T::set_gMax)
                .add_method("set_tau1", &T::set_tau1)
                .add_method("set_tau2", &T::set_tau2)
                .add_method("name", &T::name);
	    }

	    // ... add other synapse types here!
	}
    UG_REGISTRY_CATCH_THROW(grp);

}

}	// namespace ug

