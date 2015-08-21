/*
 * plugin_main.cpp
 *
 *  Created on: ??
 *      Author: pgottmann
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "lib_grid/global_attachments.h" // global attachments

#include "lib_disc/function_spaces/grid_function.h"

// Hodgkin & Huxley includes
#include "channel_interface.h"
#include "channel_hh.h"
#include "leakage.h"
#include "ElemDiscHH_base.h"
#include "ElemDiscHH_fv1.h"
#include "ElemDiscHH_Nernst_fv1.h"
#include "ElemDiscHH_Nernst_neuron_fv1.h"
#include "VM_Disc.h"
#include "cable_ass_tuner.h"
#include "order.h"

// add converted channels
#ifdef HH_CONVERTED_CHANNELS_ENABLED
#include "Convert/Debug/includefile.cpp"
#endif
//#include "Convert/Debug/hh_converted_standard_UG.h"

using namespace std;
using namespace ug::bridge;

namespace ug {
namespace cable {



/**
 *  \defgroup plugin_cable Plugin cable
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
			typedef IChannel<TDomain> T;
			string name = string("IChannel").append(suffix);
			reg.add_class_<T>(name, grp);
				//.add_method("init", &T::init)
				//.add_method("update_gating", &T::update_gating);
			reg.add_class_to_group(name, "IChannel", tag);
		}


		// Channel Interface HH
		{
			typedef ChannelHH<TDomain> T;
			typedef IChannel<TDomain> TBase;
			string name = string("ChannelHH").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductances", &T::set_conductances, "",
						"K conductance [10^6 S/m^2]| default | value=3.6e-4#"
						"Na conductance [10^6 S/m^2] | default | value=1.2e-3",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_log_mGate", &T::set_log_mGate)
				.add_method("set_log_nGate", &T::set_log_nGate)
				.add_method("set_log_hGate", &T::set_log_hGate)
				//.add_method("ionic_current", /*static_cast<void (TBase::*) (Vertex*, std::vector<double>&)> (&T::ionic_current) /*, "","", "doing flux")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelHH", tag);
		}


		// Channel Interface HH-with-Nernst
		{
			typedef ChannelHHNernst<TDomain> T;
			typedef IChannel<TDomain> TBase;
			string name = string("ChannelHHNernst").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_conductances", &T::set_conductances, "",
						"K conductance [10^6 S/m^2]| default | value=3.6e-4#"
						"Na conductance [10^6 S/m^2] | default | value=1.2e-3",
						"sets Na and K conductance values for HH mechanism")
				.add_method("set_log_mGate", &T::set_log_mGate)
				.add_method("set_log_nGate", &T::set_log_nGate)
				.add_method("set_log_hGate", &T::set_log_hGate)
				//.add_method("ionic_current", /*static_cast<void (TBase::*) (Vertex*, std::vector<double>&)> (*/&T::ionic_current) /*, "","", "doing flux")*/
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelHHNernst", tag);
		}


		// leakage "channel"
		{
			typedef ChannelLeak<TDomain> T;
			typedef IChannel<TDomain> TBase;
			string name = string("ChannelLeak").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_cond", static_cast<void (T::*)(number)>(&T::set_cond),
						"", "leak conductance in C/m^2/mV/ms| default | value=3.0e-6", "sets leak conductance for leak channel")
				.add_method("set_rev_pot", &T::set_rev_pot, "", "leakage equilibrium potential in mV | default | value=-65.0",
							"sets leakage equilibrium potential")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelLeak", tag);
		}

		// VM-Disc class
		{
			typedef IChannel<TDomain> TIChannel;
			typedef VMDisc<TDomain> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("VMDisc").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*)>
					("Subset(s)")
				.template add_constructor<void (*)(const char*, bool)>
					("Subset(s)#with ion concentrations?")
				.template add_constructor<void (*)(const char*, bool, number)>
					("Subset(s)#with ion concentrations?#InitTime")
				.add_method("set_diameter", static_cast<void (T::*)(number)>(&T::set_diameter),
						"", "new diameter | default | value=1e6", "sets a new diameter")

				.add_method("set_spec_res", static_cast<void (T::*)(number)>(&T::set_spec_res),
						"", "new specific resistance | default | value=1e6", "sets a new specific resistance")
				.add_method("set_spec_cap", static_cast<void (T::*)(number)>(&T::set_spec_cap),
						"", "new specific capacity | default | value=1e-5", "sets a new specific capacity")

				.add_method("set_k_out", static_cast<void (T::*)(number)>(&T::set_k_out),
						"", "extracellular [K] in units of mM | default | value=2.5", "sets extracellular [K]")
				.add_method("set_na_out", static_cast<void (T::*)(number)>(&T::set_na_out),
						"", "extracellular [Na] in units of mM | default | value=65.0", "sets extracellular [Na]")
				.add_method("set_ca_out", static_cast<void (T::*)(number)>(&T::set_ca_out),
						"", "extracellular [Ca] in units of mM | default | value=1.5", "sets extracellular [Ca]")

				.add_method("set_ena", static_cast<void (T::*)(number)>(&T::set_ena),
						"", "reversal potential for Na | default | value=50.0", "sets reversal potential for Na")
				.add_method("set_ek", static_cast<void (T::*)(number)>(&T::set_ek),
						"", "reversal potential for K | default | value=-77.0", "sets reversal potential for K")
				.add_method("set_eca", static_cast<void (T::*)(number)>(&T::set_eca),
						"", "reversal potential for Ca | default | value=138.0", "sets reversal potential for Ca")

				.add_method("set_temperature", static_cast<void (T::*)(number)>(&T::set_temperature),
						"", "new temperature value in units of K | default | value=310", "sets new temperature")
				.add_method("set_temperature_celsius", static_cast<void (T::*)(number)>(&T::set_temperature_celsius),
						"", "new temperature value in degrees Celsius | default | value=37", "sets new temperature")

				.add_method("set_diff_coeffs", static_cast<void (T::*)(const std::vector<number>&)> (&T::set_diff_coeffs), "",
						"diffusion coefficient of K, Na and Ca", "sets diffusion coeffizients")

				.add_method("add_channel", &T::add_channel)
				.add_method("set_influx", static_cast<void (T::*)(number, number, number, number, number, number)>(&T::set_influx), "",
						"flux value | default | value=1e-12 #"
						"x-coordinate of influx position | default | 0.0 #"
						"y-coordinate of influx position | default | 0.0 #"
						"z-coordinate of influx position | default | 0.0 #"
						"begin time | default | 0 #"
						"duration time | default | 0 ",
						"sets Position, duration, ending and influxvalue of an Influx")
				.add_method("write_gatings_for_position", &T::write_gatings_for_position)
				.add_method("set_output", static_cast<void (T::*)(bool, number, number, number, std::string)> (&T::set_output))
#ifdef UG_CPU_1
				.add_method("estimate_cfl_cond", &T::template estimate_cfl_cond<CPUAlgebra::vector_type>)
#endif
#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
				.add_method("set_synapse_handler", &T::set_synapse_handler)
#endif
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
				.add_method("set_synapse_distributor", &T::set_synapse_distributor)
#endif
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VMDisc", tag);
		}

#ifdef HH_CONVERTED_CHANNELS_ENABLED
		#include "Convert/Debug/channels.cpp"
#endif

		// Cuthill McKee ordering
		{
			reg.add_function("order_cuthillmckee", &order_cuthillmckee<TDomain>, grp.c_str(),
							 "", "approxSpace", "vertex ordering for solver optimization");
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

	}

}; // end Functionality


// end group plugin_cable
/// \}

} // end namespace cable



/// This function is called when the plugin is loaded.
extern "C" void
InitUGPlugin_HH_Kabelnew(Registry* reg, string grp)
{
	grp.append("/HH_Kabelnew");

	// declare diameter grid attachment
	typedef ANumber ADiameter;
	GlobalAttachments::declare_attachment<ADiameter>("diameter", true);

	typedef cable::Functionality Functionality;

	try
	{
		//RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}	// namespace ug

