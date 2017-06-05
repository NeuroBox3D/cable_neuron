#ifndef kdrca1_g01_converted_standard_UG_H_
#define kdrca1_g01_converted_standard_UG_H_
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
class kdrca1_g01_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
kdrca1_g01_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
        ek ( -90	*1), 
	gkdrbar(.003 *1e-05), 
        ikmax ( 0.3 *1e-05), 
        vhalfn(13   *1), 
        a0n(0.02      *1), 
        zetan(-3    *1), 
        gmn(0.7  *1), 
	nmax(2  *1), 
	q10(1*1), 
        nscale(1*1), 
m_log_SGate(false), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in kdrca1_g01_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
kdrca1_g01_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
        ek ( -90	*1), 
	gkdrbar(.003 *1e-05), 
        ikmax ( 0.3 *1e-05), 
        vhalfn(13   *1), 
        a0n(0.02      *1), 
        zetan(-3    *1), 
        gmn(0.7  *1), 
	nmax(2  *1), 
	q10(1*1), 
        nscale(1*1), 
m_log_SGate(false), 
m_log_nGate(false) {} 
UG_CATCH_THROW("Error in kdrca1_g01_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~kdrca1_g01_converted_standard_UG() {}; 
double alpn(double v); 
double betn(double v); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getek(); 
double getgkdrbar(); 
double getikmax(); 
double getvhalfn(); 
double geta0n(); 
double getzetan(); 
double getgmn(); 
double getnmax(); 
double getq10(); 
double getnscale(); 
void setek(double val); 
void setgkdrbar(double val); 
void setikmax(double val); 
void setvhalfn(double val); 
void seta0n(double val); 
void setzetan(double val); 
void setgmn(double val); 
void setnmax(double val); 
void setq10(double val); 
void setnscale(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_nGate(bool bLognGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble nGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aanGate; 
number         ek ; 
number 	gkdrbar; 
number         ikmax ; 
number         vhalfn; 
number         a0n; 
number         zetan; 
number         gmn; 
number 	nmax; 
number 	q10; 
number         nscale; 
bool m_log_SGate; 
bool m_log_nGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
const double 	gkdr = gkdrbar*n; 
const double 	ik = gkdr*(v-ek); 
}; 
 
} // namespace cable_neuron
} // namespace ug


#endif // kdrca1_g01_converted_standard_UG_H_
