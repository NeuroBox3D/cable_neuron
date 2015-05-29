#ifndef ca_converted_standard_UG_H_
#define ca_converted_standard_UG_H_
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
namespace cable {


// forward declaration 
template <typename TDomain> 
class VMDisc; 
 
template <typename TDomain> 
class ca_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
ca_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	gbar ( 0.1   	*10), 
	vshift ( 0	*1), 
	cao  ( 2.5	*1), 
cai ( 0), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1) {} 
UG_CATCH_THROW("Error in ca_converted_standard_UG initializer list. ") 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
ca_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	gbar ( 0.1   	*10), 
	vshift ( 0	*1), 
	cao  ( 2.5	*1), 
cai ( 0), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1) {} 
UG_CATCH_THROW("Error in ca_converted_standard_UG initializer list. ") 
/// destructor 
 
virtual ~ca_converted_standard_UG() {}; 
double efun(double z); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(const LocalVector& u, Edge* e); 
virtual void update_gating(number newTime, const LocalVector& u, Edge* e); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 

 
double getgbar(); 
double getvshift(); 
double getcao(); 
double getcai(); 
double gettemp(); 
double getq10(); 
double getvmin(); 
double getvmax(); 
void setgbar(double val); 
void setvshift(double val); 
void setcao(double val); 
void setcai(double val); 
void settemp(double val); 
void setq10(double val); 
void setvmin(double val); 
void setvmax(double val); 

 
protected: 
private: 
 
number m_R, m_T, m_F; 
ADouble mGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aamGate; 
ADouble hGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aahGate; 
number 	gbar ; 
number 	vshift ; 
number 	cao  ; 
number cai ; 
number 	temp ; 
number 	q10  ; 
number 	vmin ; 
number 	vmax ; 
number tadj; 
}; 
 
} // namespace cable
} // namespace ug


#endif // ca_converted_standard_UG_H_
