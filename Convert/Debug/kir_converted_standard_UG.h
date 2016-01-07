#ifndef kir_converted_standard_UG_H_
#define kir_converted_standard_UG_H_
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
class kir_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
kir_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gkbar  ( 0.00015 		*0.01), 
	mvhalf ( -52		*1), 
	mslope ( 13		*1), 
	mshift ( 30			*1), 
	qfact ( 0.5			*1), 
m_log_mGate(false) {} 
UG_CATCH_THROW("Error in kir_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
kir_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
	gkbar  ( 0.00015 		*0.01), 
	mvhalf ( -52		*1), 
	mslope ( 13		*1), 
	mshift ( 30			*1), 
	qfact ( 0.5			*1), 
m_log_mGate(false) {} 
UG_CATCH_THROW("Error in kir_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~kir_converted_standard_UG() {}; 
double taumkir(double v);
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgkbar(); 
double getmvhalf(); 
double getmslope(); 
double getmshift(); 
double getqfact(); 
void setgkbar(double val); 
void setmvhalf(double val); 
void setmslope(double val); 
void setmshift(double val); 
void setqfact(double val); 
void set_log_mGate(bool bLogmGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble mGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aamGate; 
number 	gkbar  ; 
number 	mvhalf ; 
number 	mslope ; 
number 	mshift ; 
number 	qfact ; 
bool m_log_mGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // kir_converted_standard_UG_H_
