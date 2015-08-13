#ifndef CaT_converted_standard_UG_H_
#define CaT_converted_standard_UG_H_
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
class CaT_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
CaT_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gbar ( 0.0008 *1e-05), 
	vshift ( 0	*1), 
	cao  ( 2.5	*1), 
cai ( 0), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
	v12m(50         	*1), 
	v12h(78         	*1), 
	vwm (7.4         	*1), 
	vwh(5.0         	*1), 
	am(3         	*1), 
	ah(85         	*1), 
	vm1(25         	*1), 
	vm2(100         	*1), 
	vh1(46         	*1), 
	vh2(405         	*1), 
	wm1(20         	*1), 
	wm2(15         	*1), 
	wh1(4         	*1), 
	wh2(50         	*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in CaT_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
CaT_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gbar ( 0.0008 *1e-05), 
	vshift ( 0	*1), 
	cao  ( 2.5	*1), 
cai ( 0), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
	v12m(50         	*1), 
	v12h(78         	*1), 
	vwm (7.4         	*1), 
	vwh(5.0         	*1), 
	am(3         	*1), 
	ah(85         	*1), 
	vm1(25         	*1), 
	vm2(100         	*1), 
	vh1(46         	*1), 
	vh2(405         	*1), 
	wm1(20         	*1), 
	wm2(15         	*1), 
	wh1(4         	*1), 
	wh2(50         	*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in CaT_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~CaT_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgbar(); 
double getvshift(); 
double getcao(); 
double getcai(); 
double getvmin(); 
double getvmax(); 
double getv12m(); 
double getv12h(); 
double getvwm(); 
double getvwh(); 
double getam(); 
double getah(); 
double getvm1(); 
double getvm2(); 
double getvh1(); 
double getvh2(); 
double getwm1(); 
double getwm2(); 
double getwh1(); 
double getwh2(); 
void setgbar(double val); 
void setvshift(double val); 
void setcao(double val); 
void setcai(double val); 
void setvmin(double val); 
void setvmax(double val); 
void setv12m(double val); 
void setv12h(double val); 
void setvwm(double val); 
void setvwh(double val); 
void setam(double val); 
void setah(double val); 
void setvm1(double val); 
void setvm2(double val); 
void setvh1(double val); 
void setvh2(double val); 
void setwm1(double val); 
void setwm2(double val); 
void setwh1(double val); 
void setwh2(double val); 
void set_log_mGate(bool bLogmGate); 
void set_log_hGate(bool bLoghGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble mGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aamGate; 
ADouble hGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aahGate; 
number 	gbar ; 
number 	vshift ; 
number 	cao  ; 
number cai ; 
number 	vmin ; 
number 	vmax ; 
number 	v12m; 
number 	v12h; 
number 	vwm ; 
number 	vwh; 
number 	am; 
number 	ah; 
number 	vm1; 
number 	vm2; 
number 	vh1; 
number 	vh2; 
number 	wm1; 
number 	wm2; 
number 	wh1; 
number 	wh2; 
bool m_log_mGate; 
bool m_log_hGate; 
}; 
 
} // namespace cable
} // namespace ug


#endif // CaT_converted_standard_UG_H_
