#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "bridge/util_domain_dependent.h"


#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_disc/common/local_algebra.h"
//#include "extensions/algebra_extensions.h"



// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"


// Hodgin und Huxley includes
#include "channel_interface.h"
#include "ElemDiscHH_base.h"
#include "ElemDiscHH_fv1.h"
#include "ElemDiscHH_Nernst_fv1.h"
#include "ElemDiscHH_Nernst_neuron_fv1.h"
#include "VM_Disc.h"
#include "hh_converted_standard_UG.h"

// Kabel_diff includes
//#include "kabel_diff_base.h"
//#include "kabel_diff_fe.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

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
static void Domain(bridge::Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = ug::bridge::GetDomainSuffix<TDomain>();
	string tag = ug::bridge::GetDomainTag<TDomain>();

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


	// Channel Interface hh_converted_UG
	{
		typedef hh_converted_standard_UG<TDomain> T;
		typedef IChannel<TDomain> TBase;
		string name = string("hh_UG").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)")
			//.add_method("set_conductivities", &T::set_conductivities)
			//.add_method("set_rev_pot", &T::set_rev_pot)
			//.add_method("set_accuracy", &T::set_accuracy)
			//.add_method("ionic_current", /*static_cast<void (TBase::*) (Vertex*, std::vector<double>&)> (&T::ionic_current) /*, "","", "doing flux")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "hh_UG", tag);
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
			.add_method("set_spec_res", &T::set_spec_res)
			.add_method("set_spec_cap", &T::set_spec_cap)
			.add_method("set_diff_coeffs", &T::set_diff_coeffs)
			.add_method("add_channel", &T::add_channel)
			.add_method("set_influx", &T::set_influx)
			.add_method("set_influx_ac", &T::set_influx_ac)
#ifdef PLUGIN_SYNAPSE_DISTRIBUTOR_ACTIVE
			.add_method("set_synapse_distributor", &T::set_synapse_distributor)
#endif
			//.add_method("add_func", &T::add_func)
			//.add_method("setGridFct", &T::setGridFct)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VMDisc", tag);
	}

}

}; // end Functionality



template <typename TDomain> static void
Register_Domain(bridge::Registry& reg, string grp)
{
	string suffix = bridge::GetDomainSuffix<TDomain>();
	string tag = bridge::GetDomainTag<TDomain>();


}





////////////////////////////////////////////////////////////////////////////////
//	InitUGPlugin_HH_Kabelnew
extern "C" void
InitUGPlugin_HH_Kabelnew(ug::bridge::Registry* reg, std::string parentGroup)
{
	//Registering HH-Fluxex
	std::string grpHH = parentGroup;
	grpHH.append("HH_Kabelnew/");
	//std::cout << "test: " << std::string grp(parentGroup) << endl;
	//typedef Functionality Functionality;

		//RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		//RegisterDomainDependent<Functionality>(*reg,grp);
	typedef Functionality Functionality;

	try
	{
		#ifdef UG_DIM_1
			bridge::RegisterDomain1dDependent<Functionality>(*reg, grpHH);
			//bridge::RegisterAlgebraDependent<Functionality>(*reg,grpHH);
			//bridge::RegisterDomain1dAlgebraDependent<Functionality>(*reg,grpHH);
		#endif
		#ifdef UG_DIM_2
			bridge::RegisterDomain2dDependent<Functionality>(*reg, grpHH);
			//bridge::RegisterAlgebraDependent<Functionality>(*reg,grpHH);
			//bridge::RegisterDomain2dAlgebraDependent<Functionality>(*reg,grpHH);
		#endif
		#ifdef UG_DIM_3
			bridge::RegisterDomain3dDependent<Functionality>(*reg, grpHH);
			//bridge::RegisterAlgebraDependent<Functionality>(*reg,grpHH);
			//bridge::RegisterDomain3dAlgebraDependent<Functionality>(*reg,grpHH);
		#endif
	}
	UG_REGISTRY_CATCH_THROW(grpHH);



}

}

