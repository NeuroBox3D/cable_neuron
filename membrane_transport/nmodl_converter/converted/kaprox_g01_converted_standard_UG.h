#ifndef kaprox_g01_converted_standard_UG_H_
#define kaprox_g01_converted_standard_UG_H_
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
class kaprox_g01_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
kaprox_g01_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
        ek  ( -90  *1), 
	gkabar(.008 *1e-05), 
        vhalfn(11   *1), 
        vhalfl(-56   *1), 
        a0l(0.05      *1), 
        a0n(0.05    *1), 
        zetan(-1.5    *1), 
        zetal(3    *1), 
        gmn(0.55   *1), 
        gml(1   *1), 
	lmin(2  *1), 
	nmin(0.1  *1), 
	pw(-1    *1), 
	tq(-40  *1), 
	qq(5    *1), 
	q10(5*1), 
	qtl(1*1), 
        nscale(1*1), 
        lscale(1*1), 
m_log_SGate(false), 
m_log_nGate(false), 
m_log_lGate(false) {} 
UG_CATCH_THROW("Error in kaprox_g01_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
kaprox_g01_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
        ek  ( -90  *1), 
	gkabar(.008 *1e-05), 
        vhalfn(11   *1), 
        vhalfl(-56   *1), 
        a0l(0.05      *1), 
        a0n(0.05    *1), 
        zetan(-1.5    *1), 
        zetal(3    *1), 
        gmn(0.55   *1), 
        gml(1   *1), 
	lmin(2  *1), 
	nmin(0.1  *1), 
	pw(-1    *1), 
	tq(-40  *1), 
	qq(5    *1), 
	q10(5*1), 
	qtl(1*1), 
        nscale(1*1), 
        lscale(1*1), 
m_log_SGate(false), 
m_log_nGate(false), 
m_log_lGate(false) {} 
UG_CATCH_THROW("Error in kaprox_g01_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~kaprox_g01_converted_standard_UG() {}; 
double alpn(double v); 
double betn(double v); 
double alpl(double v); 
double betl(double v); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getek(); 
double getgkabar(); 
double getvhalfn(); 
double getvhalfl(); 
double geta0l(); 
double geta0n(); 
double getzetan(); 
double getzetal(); 
double getgmn(); 
double getgml(); 
double getlmin(); 
double getnmin(); 
double getpw(); 
double gettq(); 
double getqq(); 
double getq10(); 
double getqtl(); 
double getnscale(); 
double getlscale(); 
void setek(double val); 
void setgkabar(double val); 
void setvhalfn(double val); 
void setvhalfl(double val); 
void seta0l(double val); 
void seta0n(double val); 
void setzetan(double val); 
void setzetal(double val); 
void setgmn(double val); 
void setgml(double val); 
void setlmin(double val); 
void setnmin(double val); 
void setpw(double val); 
void settq(double val); 
void setqq(double val); 
void setq10(double val); 
void setqtl(double val); 
void setnscale(double val); 
void setlscale(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_nGate(bool bLognGate); 
void set_log_lGate(bool bLoglGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble nGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aanGate; 
ADouble lGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aalGate; 
number         ek  ; 
number 	gkabar; 
number         vhalfn; 
number         vhalfl; 
number         a0l; 
number         a0n; 
number         zetan; 
number         zetal; 
number         gmn; 
number         gml; 
number 	lmin; 
number 	nmin; 
number 	pw; 
number 	tq; 
number 	qq; 
number 	q10; 
number 	qtl; 
number         nscale; 
number         lscale; 
bool m_log_SGate; 
bool m_log_nGate; 
bool m_log_lGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // kaprox_g01_converted_standard_UG_H_
