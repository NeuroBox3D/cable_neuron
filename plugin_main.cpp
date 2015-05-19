/*
 * plugin_main.cpp
 *
 *  Created on: ??
 *      Author: pgottmann
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "lib_grid/global_attachments.h" // global attachments
#include "../synapse_handler/grid/synapse_info.h"
#include "../synapse_handler/grid/synapse_info_io_traits.h"

// Hodgin und Huxley includes
#include "channel_interface.h"
#include "ElemDiscHH_base.h"
#include "ElemDiscHH_fv1.h"
#include "ElemDiscHH_Nernst_fv1.h"
#include "ElemDiscHH_Nernst_neuron_fv1.h"
#include "VM_Disc.h"
// needed to add
//#include "hh_converted_UG.h"
//#include "passive_converted_standard_UG.h"
#include "Convert/Debug/includefile.cpp"

#include <string>


using namespace std;
using namespace ug::bridge;

namespace ug {

/**
 *  This Plugin provides the Kabelgleichung in 1D
 */

/**
 * Class exporting the functionality of the plugin. All functionality that is to
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
				.add_method("set_conductivities", &T::set_conductivities)
				.add_method("set_rev_pot", &T::set_rev_pot)
				.add_method("set_accuracy", &T::set_accuracy)
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
				.add_method("set_conductivities", &T::set_conductivities)
				.add_method("set_accuracy", &T::set_accuracy)
				//.add_method("ionic_current", /*static_cast<void (TBase::*) (Vertex*, std::vector<double>&)> (*/&T::ionic_current) /*, "","", "doing flux")*/
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "ChannelHHNernst", tag);
		}


		// Only leakage Channel
		{
			typedef ChannelLeak<TDomain> T;
			typedef IChannel<TDomain> Base;
			string name = string("ChannelLeak").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
				.add_method("set_leak_cond", &T::set_leak_cond)
				.add_method("set_leak_vm", &T::set_leak_vm)
				.add_method("set_accuracy", &T::set_accuracy)
				//.add_method("ionic_current", /*static_cast<void (TBase::*) (Vertex*, std::vector<double>&)> (*/&T::ionic_current) /*, "","", "doing flux")*/
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
				.template add_constructor<void (*)(const char*,
												   SmartPtr<ApproximationSpace<TDomain> >)>
					("Subset(s)#ApproxSpace")
				.template add_constructor<void (*)(const char*,
												   SmartPtr<ApproximationSpace<TDomain> >,
												   const number)>
					("Subset(s)#ApproxSpace#InitTime")
				//.add_method("create_GridFunc", &T::create_GridFunc)
				//.add_method("getApproxSpace" , &T::getApproxSpace)
				.add_method("set_diameter", &T::set_diameter)
				//.add_method("set_diameterGeo", &T::set_diameterGeo)
				.add_method("set_spec_res", &T::set_spec_res)
				.add_method("set_spec_cap", &T::set_spec_cap)
				.add_method("set_diff_coeffs", &T::set_diff_coeffs)
				.add_method("add_channel", &T::add_channel)
				.add_method("set_influx", &T::set_influx)
				.add_method("set_ena", &T::set_ena)
				.add_method("set_ek", &T::set_ek)
				.add_method("set_eca", &T::set_eca)
				//.add_method("set_celsius", &T::set_celsius)
	#ifdef PLUGIN_SYNAPSE_HANDLER_ENABLED
				.add_method("set_synapse_handler", &T::set_synapse_handler)
	#endif
	#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED
				.add_method("set_synapse_distributor", &T::set_synapse_distributor)
	#endif
				//.add_method("add_func", &T::add_func)
				//.add_method("setGridFct", &T::setGridFct)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VMDisc", tag);
		}

//#ifdef HH_CONVERTED_CHANNELS_ENABLED
	#include "Convert/Debug/channels.cpp"
//#endif HH_CONVERTED_CHANNELS_ENABLED
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
		//static const int dim = TDomain::dim;
		string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
		string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

		typedef GridFunction<TDomain, TAlgebra> TFct;



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



////////////////////////////////////////////////////////////////////////////////
//	InitUGPlugin_HH_Kabelnew
extern "C" void
InitUGPlugin_HH_Kabelnew(Registry* reg, string grp)
{
	// declare diameter grid attachment
	typedef ANumber ADiameter;
	GlobalAttachments::declare_attachment<ADiameter>("diameter");

	//Registering HH-Fluxes
	grp.append("/HH_Kabelnew");

	try
	{
		//RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}	// namespace ug

