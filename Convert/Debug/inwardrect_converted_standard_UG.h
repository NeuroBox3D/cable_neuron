#ifndef inwardrect_converted_standard_UG_H_
#define inwardrect_converted_standard_UG_H_
#include "../../channel_interface.h" 
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
class inwardrect_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
inwardrect_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar ( 5   	*1e+06), 
	tha  ( 25	*1), 
	qa   ( 9	*1), 
	Ra   ( 0.02	*1), 
	Rb   ( 0.002	*1), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in inwardrect_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
inwardrect_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar ( 5   	*1e+06), 
	tha  ( 25	*1), 
	qa   ( 9	*1), 
	Ra   ( 0.02	*1), 
	Rb   ( 0.002	*1), 
	temp ( 23	*1), 
	q10  ( 2.3		*1), 
	vmin ( -120	*1), 
	vmax ( 100	*1), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in inwardrect_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~inwardrect_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgbar(); 
double gettha(); 
double getqa(); 
double getRa(); 
double getRb(); 
double gettemp(); 
double getq10(); 
double getvmin(); 
double getvmax(); 
void setgbar(double val); 
void settha(double val); 
void setqa(double val); 
void setRa(double val); 
void setRb(double val); 
void settemp(double val); 
void setq10(double val); 
void setvmin(double val); 
void setvmax(double val); 
void set_log_nGate(bool bLognGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble nGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aanGate; 
number 	gbar ; 
number 	tha  ; 
number 	qa   ; 
number 	Ra   ; 
number 	Rb   ; 
number 	temp ; 
number 	q10  ; 
number 	vmin ; 
number 	vmax ; 
number tadj; 
bool m_log_nGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // inwardrect_converted_standard_UG_H_
