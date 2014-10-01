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
//#include "extensions/algebra_extensions.h"



// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"


// Hodgin und Huxley includes
//#include "HH.h"
//#include "MHHDiscretization.h"
#include "ElemDiscHH_base.h"
#include "ElemDiscHH_fv1.h"

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

//	Kabel Diff Base

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
		#endif
		reg.add_class_to_group(name, "ElemDiscHH_Base", tag);
	}



//	Kabel Diff FV1
	{
		typedef ElemDiscHH_FV1<TDomain> T;
		typedef ElemDiscHH_Base<TDomain> TBase;
		string name = string("ElemDiscHH_FV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_injection", &T::set_injection)
			.add_method("set_diameter", &T::set_diameter)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ElemDiscHH_FV1", tag);
	}


}

}; // end Functionality





namespace scd {



template <typename TDomain> static void
Register_Domain(bridge::Registry& reg, string grp)
{
	string suffix = bridge::GetDomainSuffix<TDomain>();
	string tag = bridge::GetDomainTag<TDomain>();


	/*typedef HH THH;
		string name = string("HH").append(suffix);
		reg.add_class_<THH>(name, grp)
	   .add_constructor<void (*)(number, number, number)>("Spannung#Einstrom#Kapa")
	   .add_method("get_flux", &THH::get_flux,
			   	   "dt","Vm#dt","Returns HH Vm Wert.", grp)
	   .add_method("setSpannung", &THH::setSpannung,
			       "Vm#", grp)
	   .add_method("getSpannung", &THH::getSpannung,
			       "gibt spannung zurueck", grp)
	   .add_method("setEK", &THH::setEK,
			       "sets EK zurueck", grp)
	   .add_method("setEL", &THH::setEL,
				   "sets EL zurueck", grp)
	   .add_method("setENa", &THH::setENa,
				   "sets ENa zurueck", grp)
	   .add_method("setGI", &THH::setGI,
				   "sets GI zurueck", grp)
	   .add_method("setGK", &THH::setGK,
				   "sets GK zurueck", grp)
	   .add_method("setGNa", &THH::setGNa,
				   "sets GNa", grp)
	   .add_method("setEinstromPlatz", &THH::setEinstromPlatz, "SetsInjektionPlace", grp)
	   .add_method("setEinstrom", &THH::setEinstrom, "SetsInjektion", grp)
	   .add_method("setEinstromZeitDauer", &THH::setEinstromZeitDauer, "SetsInjektionDuration", grp)
	   .add_method("setEinstromWert", &THH::setEinstromWert, "SetsInjektionValue", grp)
	   .add_method("setEinstromZeit", &THH::setEinstromZeit, "SetsInjektionStarttime", grp)
	   .add_method("getEinstromPlatz", &THH::getEinstromPlatz, "GetsInjektionPlace", grp)
	   .add_method("getEinstrom", &THH::getEinstrom, "GetsInjektion", grp)
	   .add_method("getEinstromZeitDauer", &THH::getEinstromZeitDauer, "GetsInjektionDuration", grp)
	   .add_method("getEinstromWert", &THH::getEinstromWert, "GetsInjektionValue", grp)
	   .add_method("getEinstromZeit", &THH::getEinstromZeit, "GetsInjektionStarttime", grp);
	   //.set_construct_as_smart_pointer(true);
	reg.add_class_to_group(name, "HH", tag);*/


}








////////////////////////////////////////////////////////////////////////////////
//	InitUGPlugin_SynapticCalciumDynamics
extern "C" void
InitUGPlugin_HHKabelnew(ug::bridge::Registry* reg, std::string parentGroup)
{
	//Registering HH-Fluxex
	std::string grpHH = parentGroup;
	grpHH.append("HHKabelnew/");
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
		#endif
		#ifdef UG_DIM_2
			bridge::RegisterDomain2dDependent<Functionality>(*reg, grpHH);
		#endif
		#ifdef UG_DIM_3
			bridge::RegisterDomain3dDependent<Functionality>(*reg, grpHH);
		#endif
	}
	UG_REGISTRY_CATCH_THROW(grpHH);

	/*try {
		#ifdef UG_DIM_1
			bridge::RegisterDomain1dDependent<Functionality>(*reg, grpHH);
		#endif

		#ifdef UG_DIM_2
			bridge::RegisterDomain2dDependent<Functionality>(*reg, grpHH);
		#endif

		#ifdef UG_DIM_3
			bridge::RegisterDomain3dDependent<Functionality>(*reg, grpHH);
		#endif

	}*/





}

}
}

