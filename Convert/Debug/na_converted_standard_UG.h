#ifndef na_converted_standard_UG_H_
#define na_converted_standard_UG_H_
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
class na_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
na_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gbar ( 1000   	*10), 
	vshift ( -10	*1), 
	tha  ( -35	*1), 
	qa   ( 9	*1), 
	Ra   ( 0.182	*1), 
	Rb   ( 0.124	*1), 
	thi1  ( -50	*1), 
	thi2  ( -75	*1), 
	qi   ( 5	*1), 
	thinf  ( -65	*1), 
	qinf  ( 6.2	*1), 
	Rg   ( 0.0091	*1), 
	Rd   ( 0.024	*1), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in na_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
na_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gbar ( 1000   	*10), 
	vshift ( -10	*1), 
	tha  ( -35	*1), 
	qa   ( 9	*1), 
	Ra   ( 0.182	*1), 
	Rb   ( 0.124	*1), 
	thi1  ( -50	*1), 
	thi2  ( -75	*1), 
	qi   ( 5	*1), 
	thinf  ( -65	*1), 
	qinf  ( 6.2	*1), 
	Rg   ( 0.0091	*1), 
	Rd   ( 0.024	*1), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in na_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~na_converted_standard_UG() {}; 
double trap0(double v, double th, double a, double q); 
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
double gettha(); 
double getqa(); 
double getRa(); 
double getRb(); 
double getthi1(); 
double getthi2(); 
double getqi(); 
double getthinf(); 
double getqinf(); 
double getRg(); 
double getRd(); 
double gettemp(); 
double getq10(); 
double getvmin(); 
double getvmax(); 
void setgbar(double val); 
void setvshift(double val); 
void settha(double val); 
void setqa(double val); 
void setRa(double val); 
void setRb(double val); 
void setthi1(double val); 
void setthi2(double val); 
void setqi(double val); 
void setthinf(double val); 
void setqinf(double val); 
void setRg(double val); 
void setRd(double val); 
void settemp(double val); 
void setq10(double val); 
void setvmin(double val); 
void setvmax(double val); 
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
number 	tha  ; 
number 	qa   ; 
number 	Ra   ; 
number 	Rb   ; 
number 	thi1  ; 
number 	thi2  ; 
number 	qi   ; 
number 	thinf  ; 
number 	qinf  ; 
number 	Rg   ; 
number 	Rd   ; 
number 	temp ; 
number 	q10  ; 
number 	vmin ; 
number 	vmax ; 
number tadj; 
bool m_log_mGate; 
bool m_log_hGate; 
}; 
 
} // namespace cable
} // namespace ug


#endif // na_converted_standard_UG_H_
