#ifndef caL3d_converted_standard_UG_H_
#define caL3d_converted_standard_UG_H_
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
class caL3d_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
caL3d_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	p    ( 0.2e-3  	*100), 
	th   ( 5	*1), 
	q   ( 13	*1), 
	Ra   ( 1.6	*1), 
	Rb   ( 0.2	*1), 
	temp ( 22	*1), 
	q10  ( 3		*1), 
m_log_CGate(false), 
m_log_OGate(false) {} 
UG_CATCH_THROW("Error in caL3d_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
caL3d_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	p    ( 0.2e-3  	*100), 
	th   ( 5	*1), 
	q   ( 13	*1), 
	Ra   ( 1.6	*1), 
	Rb   ( 0.2	*1), 
	temp ( 22	*1), 
	q10  ( 3		*1), 
m_log_CGate(false), 
m_log_OGate(false) {} 
UG_CATCH_THROW("Error in caL3d_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~caL3d_converted_standard_UG() {}; 
double ghk(double v, double  ci, double  co); 
double efun(double z); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getp(); 
double getth(); 
double getq(); 
double getRa(); 
double getRb(); 
double gettemp(); 
double getq10(); 
void setp(double val); 
void setth(double val); 
void setq(double val); 
void setRa(double val); 
void setRb(double val); 
void settemp(double val); 
void setq10(double val); 
void set_log_CGate(bool bLogCGate); 
void set_log_OGate(bool bLogOGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble CGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaCGate; 
ADouble OGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaOGate; 
number 	p    ; 
number 	th   ; 
number 	q   ; 
number 	Ra   ; 
number 	Rb   ; 
number 	temp ; 
number 	q10  ; 
number tadj; 
bool m_log_CGate; 
bool m_log_OGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // caL3d_converted_standard_UG_H_
