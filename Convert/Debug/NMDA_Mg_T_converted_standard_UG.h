#ifndef NMDA_Mg_T_converted_standard_UG_H_
#define NMDA_Mg_T_converted_standard_UG_H_
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
class NMDA_Mg_T_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
NMDA_Mg_T_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
	Erev	( 5    	*1), 
	gmax	( 500  	*1e+15), 
	mg	( 1  	*1), 
	vmin 	( -120	*1), 
	vmax 	( 100	*1), 
	valence ( -2	*1), 
	memb_fraction ( 0.8*1), 
	Rb		( 10e-3    	*1), 
	Ru		( 0.02016 	*1), 
	Ro		( 46.5e-3   	*1), 
	Rc		( 91.6e-3   	*1), 
	Rd1		( 0.02266  	*1), 
	Rr1		( 0.00736  	*1), 
	Rd2 		( 0.004429	*1), 
	Rr2 		( 0.0023	*1), 
	Rmb		( 0.05e-3	*1), 
	Rmu		( 12800e-3	*1), 
	Rmc1b		( 0.00005e-3	*1), 
	Rmc1u		( 0.06	*1), 
	Rmc2b		( 0.00005e-3	*1), 
	Rmc2u		( 0.06	*1), 
	Rmd1b		( 0.00005e-3	*1), 
	Rmd1u		( 0.06	*1), 
	Rmd2b		( 0.00005e-3	*1), 
	Rmd2u		( 0.06	*1), 
	RbMg		( 10e-3		*1), 
	RuMg		( 0.06156	*1), 
	RoMg		( 46.5e-3		*1), 
	RcMg		( 91.6e-3	*1), 
	Rd1Mg		( 0.02163	*1), 
	Rr1Mg		( 0.004002	*1), 
	Rd2Mg		( 0.002678	*1), 
	Rr2Mg		( 0.001932	*1), 
C (0), 
m_log_UGate(false), 
m_log_FGate(false), 
m_log_ClGate(false), 
m_log_D1Gate(false), 
m_log_D2Gate(false), 
m_log_OGate(false), 
m_log_UMgGate(false), 
m_log_ClMgGate(false), 
m_log_D1MgGate(false), 
m_log_D2MgGate(false), 
m_log_OMgGate(false) {} 
UG_CATCH_THROW("Error in NMDA_Mg_T_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
NMDA_Mg_T_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
	Erev	( 5    	*1), 
	gmax	( 500  	*1e+15), 
	mg	( 1  	*1), 
	vmin 	( -120	*1), 
	vmax 	( 100	*1), 
	valence ( -2	*1), 
	memb_fraction ( 0.8*1), 
	Rb		( 10e-3    	*1), 
	Ru		( 0.02016 	*1), 
	Ro		( 46.5e-3   	*1), 
	Rc		( 91.6e-3   	*1), 
	Rd1		( 0.02266  	*1), 
	Rr1		( 0.00736  	*1), 
	Rd2 		( 0.004429	*1), 
	Rr2 		( 0.0023	*1), 
	Rmb		( 0.05e-3	*1), 
	Rmu		( 12800e-3	*1), 
	Rmc1b		( 0.00005e-3	*1), 
	Rmc1u		( 0.06	*1), 
	Rmc2b		( 0.00005e-3	*1), 
	Rmc2u		( 0.06	*1), 
	Rmd1b		( 0.00005e-3	*1), 
	Rmd1u		( 0.06	*1), 
	Rmd2b		( 0.00005e-3	*1), 
	Rmd2u		( 0.06	*1), 
	RbMg		( 10e-3		*1), 
	RuMg		( 0.06156	*1), 
	RoMg		( 46.5e-3		*1), 
	RcMg		( 91.6e-3	*1), 
	Rd1Mg		( 0.02163	*1), 
	Rr1Mg		( 0.004002	*1), 
	Rd2Mg		( 0.002678	*1), 
	Rr2Mg		( 0.001932	*1), 
C (0), 
m_log_UGate(false), 
m_log_FGate(false), 
m_log_ClGate(false), 
m_log_D1Gate(false), 
m_log_D2Gate(false), 
m_log_OGate(false), 
m_log_UMgGate(false), 
m_log_ClMgGate(false), 
m_log_D1MgGate(false), 
m_log_D2MgGate(false), 
m_log_OMgGate(false) {} 
UG_CATCH_THROW("Error in NMDA_Mg_T_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~NMDA_Mg_T_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void approx_space_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getErev(); 
double getgmax(); 
double getmg(); 
double getvmin(); 
double getvmax(); 
double getvalence(); 
double getmemb_fraction(); 
double getRb(); 
double getRu(); 
double getRo(); 
double getRc(); 
double getRd1(); 
double getRr1(); 
double getRd2(); 
double getRr2(); 
double getRmb(); 
double getRmu(); 
double getRmc1b(); 
double getRmc1u(); 
double getRmc2b(); 
double getRmc2u(); 
double getRmd1b(); 
double getRmd1u(); 
double getRmd2b(); 
double getRmd2u(); 
double getRbMg(); 
double getRuMg(); 
double getRoMg(); 
double getRcMg(); 
double getRd1Mg(); 
double getRr1Mg(); 
double getRd2Mg(); 
double getRr2Mg(); 
double getC(); 
void setErev(double val); 
void setgmax(double val); 
void setmg(double val); 
void setvmin(double val); 
void setvmax(double val); 
void setvalence(double val); 
void setmemb_fraction(double val); 
void setRb(double val); 
void setRu(double val); 
void setRo(double val); 
void setRc(double val); 
void setRd1(double val); 
void setRr1(double val); 
void setRd2(double val); 
void setRr2(double val); 
void setRmb(double val); 
void setRmu(double val); 
void setRmc1b(double val); 
void setRmc1u(double val); 
void setRmc2b(double val); 
void setRmc2u(double val); 
void setRmd1b(double val); 
void setRmd1u(double val); 
void setRmd2b(double val); 
void setRmd2u(double val); 
void setRbMg(double val); 
void setRuMg(double val); 
void setRoMg(double val); 
void setRcMg(double val); 
void setRd1Mg(double val); 
void setRr1Mg(double val); 
void setRd2Mg(double val); 
void setRr2Mg(double val); 
void setC(double val); 
void set_log_UGate(bool bLogUGate); 
void set_log_FGate(bool bLogFGate); 
void set_log_ClGate(bool bLogClGate); 
void set_log_D1Gate(bool bLogD1Gate); 
void set_log_D2Gate(bool bLogD2Gate); 
void set_log_OGate(bool bLogOGate); 
void set_log_UMgGate(bool bLogUMgGate); 
void set_log_ClMgGate(bool bLogClMgGate); 
void set_log_D1MgGate(bool bLogD1MgGate); 
void set_log_D2MgGate(bool bLogD2MgGate); 
void set_log_OMgGate(bool bLogOMgGate); 

 
protected: 
private: 
 
ADouble UGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaUGate; 
ADouble FGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaFGate; 
ADouble ClGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaClGate; 
ADouble D1Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaD1Gate; 
ADouble D2Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaD2Gate; 
ADouble OGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaOGate; 
ADouble UMgGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaUMgGate; 
ADouble ClMgGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaClMgGate; 
ADouble D1MgGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaD1MgGate; 
ADouble D2MgGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaD2MgGate; 
ADouble OMgGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaOMgGate; 
number 	Erev	; 
number 	gmax	; 
number 	mg	; 
number 	vmin 	; 
number 	vmax 	; 
number 	valence ; 
number 	memb_fraction ; 
number 	Rb		; 
number 	Ru		; 
number 	Ro		; 
number 	Rc		; 
number 	Rd1		; 
number 	Rr1		; 
number 	Rd2 		; 
number 	Rr2 		; 
number 	Rmb		; 
number 	Rmu		; 
number 	Rmc1b		; 
number 	Rmc1u		; 
number 	Rmc2b		; 
number 	Rmc2u		; 
number 	Rmd1b		; 
number 	Rmd1u		; 
number 	Rmd2b		; 
number 	Rmd2u		; 
number 	RbMg		; 
number 	RuMg		; 
number 	RoMg		; 
number 	RcMg		; 
number 	Rd1Mg		; 
number 	Rr1Mg		; 
number 	Rd2Mg		; 
number 	Rr2Mg		; 
number C ; 
bool m_log_UGate; 
bool m_log_FGate; 
bool m_log_ClGate; 
bool m_log_D1Gate; 
bool m_log_D2Gate; 
bool m_log_OGate; 
bool m_log_UMgGate; 
bool m_log_ClMgGate; 
bool m_log_D1MgGate; 
bool m_log_D2MgGate; 
bool m_log_OMgGate; 
}; 
 
} // namespace cable
} // namespace ug


#endif // NMDA_Mg_T_converted_standard_UG_H_
