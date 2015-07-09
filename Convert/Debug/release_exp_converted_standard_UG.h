#ifndef release_exp_converted_standard_UG_H_
#define release_exp_converted_standard_UG_H_
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
class release_exp_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
release_exp_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
    tau1(.1 *1), 
    tau2 ( 10 *1), 
m_log_AGate(false), 
m_log_BGate(false) {} 
UG_CATCH_THROW("Error in release_exp_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
release_exp_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
    tau1(.1 *1), 
    tau2 ( 10 *1), 
m_log_AGate(false), 
m_log_BGate(false) {} 
UG_CATCH_THROW("Error in release_exp_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~release_exp_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double gettau1(); 
double gettau2(); 
void settau1(double val); 
void settau2(double val); 
void set_log_AGate(bool bLogAGate); 
void set_log_BGate(bool bLogBGate); 

 
protected: 
private: 
 
ADouble AGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaAGate; 
ADouble BGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaBGate; 
number     tau1; 
number     tau2 ; 
number total; 
bool m_log_AGate; 
bool m_log_BGate; 
}; 
 
} // namespace cable
} // namespace ug


#endif // release_exp_converted_standard_UG_H_
