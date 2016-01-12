#ifndef release_exp_converted_standard_UG_H_
#define release_exp_converted_standard_UG_H_
#include "../../cable_membrane_transport_interface.h" 
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

#include "../../../cable_disc/cable_equation.h" 
 
#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class release_exp_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
release_exp_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
    tau1(.1 *1), 
    tau2 ( 10 *1), 
m_log_SGate(false), 
m_log_AGate(false), 
m_log_BGate(false) {} 
UG_CATCH_THROW("Error in release_exp_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
release_exp_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
    tau1(.1 *1), 
    tau2 ( 10 *1), 
m_log_SGate(false), 
m_log_AGate(false), 
m_log_BGate(false) {} 
UG_CATCH_THROW("Error in release_exp_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~release_exp_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double gettau1(); 
double gettau2(); 
void settau1(double val); 
void settau2(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_AGate(bool bLogAGate); 
void set_log_BGate(bool bLogBGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble AGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaAGate; 
ADouble BGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaBGate; 
number     tau1; 
number     tau2 ; 
number total; 
bool m_log_SGate; 
bool m_log_AGate; 
bool m_log_BGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
const double     tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1); 
const double     factor = 1/(-exp(-tp/tau1) + exp(-tp/tau2));

}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // release_exp_converted_standard_UG_H_