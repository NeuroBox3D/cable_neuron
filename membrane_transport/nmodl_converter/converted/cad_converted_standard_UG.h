#ifndef cad_converted_standard_UG_H_
#define cad_converted_standard_UG_H_
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

#include "../../cable_membrane_transport_interface.h"
#include "bindings/lua/lua_user_data.h" 
namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class cad_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
cad_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	depth	( .1	*1e+06), 
	taur	( 200	*1), 
	cainf	( 100e-6*1), 
cai ( 0), 
m_log_SGate(false), 
m_log_caSGate(false) {} 
UG_CATCH_THROW("Error in cad_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
cad_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	depth	( .1	*1e+06), 
	taur	( 200	*1), 
	cainf	( 100e-6*1), 
cai ( 0), 
m_log_SGate(false), 
m_log_caSGate(false) {} 
UG_CATCH_THROW("Error in cad_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~cad_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getdepth(); 
double gettaur(); 
double getcainf(); 
double getcai(); 
void setdepth(double val); 
void settaur(double val); 
void setcainf(double val); 
void setcai(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_caSGate(bool bLogcaSGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble caSGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aacaSGate; 
number 	depth	; 
number 	taur	; 
number 	cainf	; 
number cai ; 
bool m_log_SGate; 
bool m_log_caSGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
const double 	ca = cainf; 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // cad_converted_standard_UG_H_
