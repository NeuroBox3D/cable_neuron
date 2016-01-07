#ifndef ar_converted_standard_UG_H_
#define ar_converted_standard_UG_H_
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

#include "../../cable_membrane_transport_interface.h"
#include "bindings/lua/lua_user_data.h" 
namespace ug {
namespace cable {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class ar_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
ar_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
         g0 ( .0001      *0.01), 
         e ( -75 	  	 *1), 
         c ( 1000000 	 *1e+11){} 
UG_CATCH_THROW("Error in ar_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
ar_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
         g0 ( .0001      *0.01), 
         e ( -75 	  	 *1), 
         c ( 1000000 	 *1e+11){} 
UG_CATCH_THROW("Error in ar_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~ar_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getg0(); 
double gete(); 
double getc(); 
void setg0(double val); 
void sete(double val); 
void setc(double val); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
number          g0 ; 
number          e ; 
number          c ; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // ar_converted_standard_UG_H_
