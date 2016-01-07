#ifndef h_converted_standard_UG_H_
#define h_converted_standard_UG_H_
#include "../../channel_interface.h" 
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
namespace cable {


// forward declaration 
template <typename TDomain> 
class CableEquation; 
 
template <typename TDomain> 
class h_converted_standard_UG
    : public ICableMembraneTransport<TDomain> 
{ 
    public: 
 using ICableMembraneTransport <TDomain>::m_pCE; 
 

 
 



/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(cont char*) 
h_converted_standard_UG(const char* functions, const char* subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
ehd ( 0), 
    ghdbar(.0001 *1e-05), 
    vhalfl(-81   *1), 
    kl(-8*1), 
    vhalft(-75   *1), 
    a0t(0.011      *1), 
    zetat(2.2    *1), 
    gmt(.4   *1), 
    q10(4.5*1), 
    qtl(1*1), 
m_log_SGate(false), 
m_log_lGate(false) {} 
UG_CATCH_THROW("Error in h_converted_standard_UG initializer list. "); 
 
 
/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&) 
h_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : ICableMembraneTransport<TDomain>(functions, subsets), 
ehd ( 0), 
    ghdbar(.0001 *1e-05), 
    vhalfl(-81   *1), 
    kl(-8*1), 
    vhalft(-75   *1), 
    a0t(0.011      *1), 
    zetat(2.2    *1), 
    gmt(.4   *1), 
    q10(4.5*1), 
    qtl(1*1), 
m_log_SGate(false), 
m_log_lGate(false) {} 
UG_CATCH_THROW("Error in h_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~h_converted_standard_UG() {}; 
double alpt(double v); 
double bett(double v); 
/// create attachments and accessors 
void init_attachments(); 
// inherited from ICableMembraneTransport 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void ce_obj_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getehd(); 
double getghdbar(); 
double getvhalfl(); 
double getkl(); 
double getvhalft(); 
double geta0t(); 
double getzetat(); 
double getgmt(); 
double getq10(); 
double getqtl(); 
void setehd(double val); 
void setghdbar(double val); 
void setvhalfl(double val); 
void setkl(double val); 
void setvhalft(double val); 
void seta0t(double val); 
void setzetat(double val); 
void setgmt(double val); 
void setq10(double val); 
void setqtl(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_lGate(bool bLoglGate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble lGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aalGate; 
number ehd ; 
number     ghdbar; 
number     vhalfl; 
number     kl; 
number     vhalft; 
number     a0t; 
number     zetat; 
number     gmt; 
number     q10; 
number     qtl; 
number taul; 
number linf; 
bool m_log_SGate; 
bool m_log_lGate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // h_converted_standard_UG_H_
