#ifndef nax_g01_converted_standard_UG_H_
#define nax_g01_converted_standard_UG_H_
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
class nax_g01_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
nax_g01_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar ( 0.010   	*1e-05), 
	tha  (  -30	*1), 
	qa   ( 7.2	*1), 
	Ra   ( 0.4	*1), 
	Rb   ( 0.124 	*1), 
	thi1  ( -45	*1), 
	thi2  ( -45 	*1), 
	qd   ( 1.5	*1), 
	qg   ( 1.5      *1), 
	mmin(0.02	*1), 
	hmin(0.5	*1), 
	q10(2*1), 
	Rg   ( 0.01 	*1), 
	Rd   ( .03 	*1), 
	thinf  ( -50 	*1), 
	qinf  ( 4 	*1), 
	ena ( 55	*1), 
        mscale ( 1*1), 
        hscale ( 1*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in nax_g01_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
nax_g01_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
	gbar ( 0.010   	*1e-05), 
	tha  (  -30	*1), 
	qa   ( 7.2	*1), 
	Ra   ( 0.4	*1), 
	Rb   ( 0.124 	*1), 
	thi1  ( -45	*1), 
	thi2  ( -45 	*1), 
	qd   ( 1.5	*1), 
	qg   ( 1.5      *1), 
	mmin(0.02	*1), 
	hmin(0.5	*1), 
	q10(2*1), 
	Rg   ( 0.01 	*1), 
	Rd   ( .03 	*1), 
	thinf  ( -50 	*1), 
	qinf  ( 4 	*1), 
	ena ( 55	*1), 
        mscale ( 1*1), 
        hscale ( 1*1), 
m_log_mGate(false), 
m_log_hGate(false) {} 
UG_CATCH_THROW("Error in nax_g01_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~nax_g01_converted_standard_UG() {}; 
double trap0(double v, double th, double a, double q); 
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
double getthi1(); 
double getthi2(); 
double getqd(); 
double getqg(); 
double getmmin(); 
double gethmin(); 
double getq10(); 
double getRg(); 
double getRd(); 
double getthinf(); 
double getqinf(); 
double getena(); 
double getmscale(); 
double gethscale(); 
void setgbar(double val); 
void settha(double val); 
void setqa(double val); 
void setRa(double val); 
void setRb(double val); 
void setthi1(double val); 
void setthi2(double val); 
void setqd(double val); 
void setqg(double val); 
void setmmin(double val); 
void sethmin(double val); 
void setq10(double val); 
void setRg(double val); 
void setRd(double val); 
void setthinf(double val); 
void setqinf(double val); 
void setena(double val); 
void setmscale(double val); 
void sethscale(double val); 
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
number 	tha  ; 
number 	qa   ; 
number 	Ra   ; 
number 	Rb   ; 
number 	thi1  ; 
number 	thi2  ; 
number 	qd   ; 
number 	qg   ; 
number 	mmin; 
number 	hmin; 
number 	q10; 
number 	Rg   ; 
number 	Rd   ; 
number 	thinf  ; 
number 	qinf  ; 
number 	ena ; 
number         mscale ; 
number         hscale ; 
bool m_log_mGate; 
bool m_log_hGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
const double         thegna = gbar*m*m*m*h; 
const double 	ina = thegna * (v - ena); 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // nax_g01_converted_standard_UG_H_
