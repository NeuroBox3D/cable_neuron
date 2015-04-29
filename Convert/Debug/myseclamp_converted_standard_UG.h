#ifndef myseclamp_converted_standard_UG_H_
#define myseclamp_converted_standard_UG_H_
#include "channel_interface.h" 
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

#include "../../VM_Disc.h" 
 
#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug { 
 
// forward declaration 
template <typename TDomain> 
class VMDisc; 
 
template <typename TDomain> 
class myseclamp_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
myseclamp_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	rs ( 1 *1), 
dur1 ( 0), 
dur2 ( 0), 
dur3 ( 0), 
dur4 ( 0) {} 
UG_CATCH_THROW("Error in myseclamp_converted_standard_UG initializer list. ") 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
myseclamp_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	rs ( 1 *1), 
dur1 ( 0), 
dur2 ( 0), 
dur3 ( 0), 
dur4 ( 0) {} 
UG_CATCH_THROW("Error in myseclamp_converted_standard_UG initializer list. ") 
/// destructor 
 
virtual ~myseclamp_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(const LocalVector& u, Edge* e); 
virtual void update_gating(number newTime, const LocalVector& u, Edge* e); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 

 
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
 
number m_R, m_T, m_F; 
number 	rs ; 
number dur1 ; 
number dur2 ; 
number dur3 ; 
number dur4 ; 
}; 
 
} // namespace ug 
 
 
#endif // myseclamp_converted_standard_UG_H_
