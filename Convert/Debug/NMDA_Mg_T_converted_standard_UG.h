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

#include "VM_Disc.h" 
 
#include <vector> 
#include <stdio.h> 
#include "bindings/lua/lua_user_data.h" 
namespace ug { 
 
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
m_R(8.314), m_T(293.0), m_F(96485.0), 
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
	Rr2Mg		( 0.001932	*1) {} 
UG_CATCH_THROW("Error in NMDA_Mg_T_converted_standard_UG initializer list. ") 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
NMDA_Mg_T_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
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
	Rr2Mg		( 0.001932	*1) {} 
UG_CATCH_THROW("Error in NMDA_Mg_T_converted_standard_UG initializer list. ") 
/// destructor 
 
virtual ~NMDA_Mg_T_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(const LocalVector& u, Edge* e); 
virtual void update_gating(number newTime, const LocalVector& u, Edge* e); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 

 
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
number memb_fraction; 
number vmin; 
number  vmax; 
number  valence; 
number 	Erev	; 
number 	gmax	; 
number 	mg	; 
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
number m_R, m_T, m_F; 
}; 
 
} // namespace ug 
 
 
#endif // NMDA_Mg_T_converted_standard_UG_H_
