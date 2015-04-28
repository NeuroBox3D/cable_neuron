#ifndef caL3d_converted_standard_UG_H_
#define caL3d_converted_standard_UG_H_
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
class caL3d_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
caL3d_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	p    ( 0.2e-3  	*100), 
	th   ( 5	*1), 
	q   ( 13	*1), 
	Ra   ( 1.6	*1), 
	Rb   ( 0.2	*1), 
	temp ( 22	*1), 
	q10  ( 3		*1) {} 
UG_CATCH_THROW("Error in caL3d_converted_standard_UG initializer list. ") 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
caL3d_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
	p    ( 0.2e-3  	*100), 
	th   ( 5	*1), 
	q   ( 13	*1), 
	Ra   ( 1.6	*1), 
	Rb   ( 0.2	*1), 
	temp ( 22	*1), 
	q10  ( 3		*1) {} 
UG_CATCH_THROW("Error in caL3d_converted_standard_UG initializer list. ") 
/// destructor 
 
virtual ~caL3d_converted_standard_UG() {}; 
double ghk(double v, double  ci, double  co); 
double efun(double z); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(const LocalVector& u, Edge* e); 
virtual void update_gating(number newTime, const LocalVector& u, Edge* e); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 

 
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

 
protected: 
private: 
 
number m_R, m_T, m_F; 
ADouble CGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaCGate; 
ADouble OGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaOGate; 
number tadj; 
number q10; 
number temp; 
number 	p    ; 
number 	th   ; 
number 	q   ; 
number 	Ra   ; 
number 	Rb   ; 
}; 
 
} // namespace ug 
 
 
#endif // caL3d_converted_standard_UG_H_
