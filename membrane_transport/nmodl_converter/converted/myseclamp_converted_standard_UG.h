#ifndef myseclamp_converted_standard_UG_H_
#define myseclamp_converted_standard_UG_H_
#include "../../cable_membrane_transport_interface.h" 
#include "lib_grid/lg_base.h" 
#include "lib_grid/grid/grid_base_objects.h" 

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h" 
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h" 
#include "lib_disc/spatial_disc/disc_util/geom_provider.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include "lib_disc/common/local_algebra.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/interpolate.h" 

#include "bridge/bridge.h" 
#include "bridge/util.h" 
#include "bridge/util_domain_algebra_dependent.h" 
#include "bridge/util_domain_dependent.h" 

#include "common/util/smart_pointer.h" 
#include "common/util/vector_util.h" 

#include "../../../cable_disc/cable_equation.h" 
 
#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class myseclamp_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
myseclamp_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	rs ( 1 *1), 
dur1 ( 0), 
dur2 ( 0), 
dur3 ( 0), 
dur4 ( 0){} 
UG_CATCH_THROW("Error in myseclamp_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
myseclamp_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	rs ( 1 *1), 
dur1 ( 0), 
dur2 ( 0), 
dur3 ( 0), 
dur4 ( 0){} 
UG_CATCH_THROW("Error in myseclamp_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~myseclamp_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getrs(); 
double getdur1(); 
double getdur2(); 
double getdur3(); 
double getdur4(); 
void setrs(double val); 
void setdur1(double val); 
void setdur2(double val); 
void setdur3(double val); 
void setdur4(double val); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
number 	rs ; 
number dur1 ; 
number dur2 ; 
number dur3 ; 
number dur4 ; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
const double 	tc2 = dur1 + dur2; 
const double 	tc3 = tc2 + dur3; 
const double 	tc4 = tc3 + dur4; 
const double 	on = 0; 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // myseclamp_converted_standard_UG_H_