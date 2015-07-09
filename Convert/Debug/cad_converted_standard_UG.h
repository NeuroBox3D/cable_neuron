#ifndef cad_converted_standard_UG_H_
#define cad_converted_standard_UG_H_
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
class cad_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
cad_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
	depth	( .1	*1e+06), 
	taur	( 200	*1), 
	cainf	( 100e-6*1), 
cai ( 0), 
m_log_caSGate(false) {} 
UG_CATCH_THROW("Error in cad_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
cad_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
	depth	( .1	*1e+06), 
	taur	( 200	*1), 
	cainf	( 100e-6*1), 
cai ( 0), 
m_log_caSGate(false) {} 
UG_CATCH_THROW("Error in cad_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~cad_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getdepth(); 
double gettaur(); 
double getcainf(); 
double getcai(); 
void setdepth(double val); 
void settaur(double val); 
void setcainf(double val); 
void setcai(double val); 
void set_log_caSGate(bool bLogcaSGate); 

 
protected: 
private: 
 
ADouble caSGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aacaSGate; 
number 	depth	; 
number 	taur	; 
number 	cainf	; 
number cai ; 
bool m_log_caSGate; 
}; 
 
} // namespace cable
} // namespace ug


#endif // cad_converted_standard_UG_H_
