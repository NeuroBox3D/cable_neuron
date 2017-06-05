#ifndef h_g05_converted_standard_UG_H_
#define h_g05_converted_standard_UG_H_
#include "../../membrane_transport/cable_membrane_transport_interface.h" 
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
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class h_g05_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
h_g05_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar(8.0			*1e+06), 
	erevh(-13			*1), 
	vhalf(-81			*1), 
    a0(0.00057			*1), 
	zeta(7				*1), 
    ab(0.4   			*1), 
	qten(4.5			*1), 
    temp ( 33			*1), 
	gas ( 8.315			*1), 
	farad ( 9.648e4		*1), 
m_log_SGate(false), 
m_log_hhGate(false) {} 
UG_CATCH_THROW("Error in h_g05_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
h_g05_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar(8.0			*1e+06), 
	erevh(-13			*1), 
	vhalf(-81			*1), 
    a0(0.00057			*1), 
	zeta(7				*1), 
    ab(0.4   			*1), 
	qten(4.5			*1), 
    temp ( 33			*1), 
	gas ( 8.315			*1), 
	farad ( 9.648e4		*1), 
m_log_SGate(false), 
m_log_hhGate(false) {} 
UG_CATCH_THROW("Error in h_g05_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~h_g05_converted_standard_UG() {}; 
double alpha(double v); 
double beta(double v); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgbar(); 
double geterevh(); 
double getvhalf(); 
double geta0(); 
double getzeta(); 
double getab(); 
double getqten(); 
double gettemp(); 
double getgas(); 
double getfarad(); 
void setgbar(double val); 
void seterevh(double val); 
void setvhalf(double val); 
void seta0(double val); 
void setzeta(double val); 
void setab(double val); 
void setqten(double val); 
void settemp(double val); 
void setgas(double val); 
void setfarad(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_hhGate(bool bLoghhGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble hhGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aahhGate; 
number 	gbar; 
number 	erevh; 
number 	vhalf; 
number     a0; 
number 	zeta; 
number     ab; 
number 	qten; 
number     temp ; 
number 	gas ; 
number 	farad ; 
number tau; 
number inf; 
bool m_log_SGate; 
bool m_log_hhGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // h_g05_converted_standard_UG_H_
