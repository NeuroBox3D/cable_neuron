#ifndef hh_converted_standard_UG_H_
#define hh_converted_standard_UG_H_
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
class hh_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
hh_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
        gnabar ( .12 *0.01), 
        gkbar ( .036 *0.01), 
        gl ( .0003 *0.01), 
        el ( -54.3 *1), 
m_log_SGate(false), 
m_log_mGate(false), 
m_log_hGate(false), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in hh_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
hh_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
        gnabar ( .12 *0.01), 
        gkbar ( .036 *0.01), 
        gl ( .0003 *0.01), 
        el ( -54.3 *1), 
m_log_SGate(false), 
m_log_mGate(false), 
m_log_hGate(false), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in hh_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~hh_converted_standard_UG() {}; 
double vtrap(double x, double y); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgnabar(); 
double getgkbar(); 
double getgl(); 
double getel(); 
void setgnabar(double val); 
void setgkbar(double val); 
void setgl(double val); 
void setel(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_mGate(bool bLogmGate); 
void set_log_hGate(bool bLoghGate); 
void set_log_nGate(bool bLognGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble mGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aamGate; 
ADouble hGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aahGate; 
ADouble nGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aanGate; 
number         gnabar ; 
number         gkbar ; 
number         gl ; 
number         el ; 
number ntau; 
number minf; 
number hinf; 
number ninf; 
number mtau; 
number htau; 
bool m_log_SGate; 
bool m_log_mGate; 
bool m_log_hGate; 
bool m_log_nGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // hh_converted_standard_UG_H_
