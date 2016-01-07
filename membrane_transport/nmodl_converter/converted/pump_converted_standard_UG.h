#ifndef pump_converted_standard_UG_H_
#define pump_converted_standard_UG_H_
#include "../../membrane_transport/cable_membrane_transport_interface.h" 
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

#include "../../cable_equation.h" 
 
#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug {
namespace cable {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class pump_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
pump_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
nai ( 0), 
        ipumpmax  ( 0.0036   *1e-05), 
        km ( 10.0        *1), 
        n(1.5*1), 
        nainit ( 4  *1){} 
UG_CATCH_THROW("Error in pump_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
pump_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
nai ( 0), 
        ipumpmax  ( 0.0036   *1e-05), 
        km ( 10.0        *1), 
        n(1.5*1), 
        nainit ( 4  *1){} 
UG_CATCH_THROW("Error in pump_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~pump_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getnai(); 
double getipumpmax(); 
double getkm(); 
double getn(); 
double getnainit(); 
void setnai(double val); 
void setipumpmax(double val); 
void setkm(double val); 
void setn(double val); 
void setnainit(double val); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
number nai ; 
number         ipumpmax  ; 
number         km ; 
number         n; 
number         nainit ; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // pump_converted_standard_UG_H_
