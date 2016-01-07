#ifndef Kv4_csiosi_converted_standard_UG_H_
#define Kv4_csiosi_converted_standard_UG_H_
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
class Kv4_csiosi_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
Kv4_csiosi_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
      	gmax ( 0.00			*0.01), 
      	ek(-90              *1), 
		a ( 7				*1), 
		za ( 0.315646648				*1), 
      	b ( .090			*1), 
      	zb ( -2.062276				*1), 
		c ( 1.01216107		*0.001), 
		zc ( 0.500095665*1), 
		d ( 2.498881		*1), 
		zd ( -1.1546687*1), 
		k ( 7.69049072		*1), 
		zk ( 0.05502051*1), 
		l ( 4.38562354		*1), 
		zl ( -0.07092366 *1), 
		f ( 0.277130485			*1), 
		q ( 1.01314807			*1), 
		kci ( 0.121900093	*1), 
		kic ( 0.0017935468 	*1), 
		koi ( 0.5195		*1), 
		kio ( 0.0436		*1), 
		kii2 ( 0.1475		*1), 
		ki2i ( 0.0333		*1), 
m_log_SGate(false), 
m_log_C0Gate(false), 
m_log_C1Gate(false), 
m_log_C2Gate(false), 
m_log_C3Gate(false), 
m_log_C4Gate(false), 
m_log_C5Gate(false), 
m_log_I0Gate(false), 
m_log_I1Gate(false), 
m_log_I2Gate(false), 
m_log_I3Gate(false), 
m_log_I4Gate(false), 
m_log_I5Gate(false), 
m_log_OGate(false), 
m_log_I6Gate(false), 
m_log_I7Gate(false) {} 
UG_CATCH_THROW("Error in Kv4_csiosi_converted_standard_UG initializer list. "); 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
Kv4_csiosi_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
      	gmax ( 0.00			*0.01), 
      	ek(-90              *1), 
		a ( 7				*1), 
		za ( 0.315646648				*1), 
      	b ( .090			*1), 
      	zb ( -2.062276				*1), 
		c ( 1.01216107		*0.001), 
		zc ( 0.500095665*1), 
		d ( 2.498881		*1), 
		zd ( -1.1546687*1), 
		k ( 7.69049072		*1), 
		zk ( 0.05502051*1), 
		l ( 4.38562354		*1), 
		zl ( -0.07092366 *1), 
		f ( 0.277130485			*1), 
		q ( 1.01314807			*1), 
		kci ( 0.121900093	*1), 
		kic ( 0.0017935468 	*1), 
		koi ( 0.5195		*1), 
		kio ( 0.0436		*1), 
		kii2 ( 0.1475		*1), 
		ki2i ( 0.0333		*1), 
m_log_SGate(false), 
m_log_C0Gate(false), 
m_log_C1Gate(false), 
m_log_C2Gate(false), 
m_log_C3Gate(false), 
m_log_C4Gate(false), 
m_log_C5Gate(false), 
m_log_I0Gate(false), 
m_log_I1Gate(false), 
m_log_I2Gate(false), 
m_log_I3Gate(false), 
m_log_I4Gate(false), 
m_log_I5Gate(false), 
m_log_OGate(false), 
m_log_I6Gate(false), 
m_log_I7Gate(false) {} 
UG_CATCH_THROW("Error in Kv4_csiosi_converted_standard_UG initializer list. "); 
/// destructor 
 
virtual ~Kv4_csiosi_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 
virtual std::vector<number> state_values(number x, number y, number z); 

 
double getgmax(); 
double getek(); 
double geta(); 
double getza(); 
double getb(); 
double getzb(); 
double getc(); 
double getzc(); 
double getd(); 
double getzd(); 
double getk(); 
double getzk(); 
double getl(); 
double getzl(); 
double getf(); 
double getq(); 
double getkci(); 
double getkic(); 
double getkoi(); 
double getkio(); 
double getkii2(); 
double getki2i(); 
void setgmax(double val); 
void setek(double val); 
void seta(double val); 
void setza(double val); 
void setb(double val); 
void setzb(double val); 
void setc(double val); 
void setzc(double val); 
void setd(double val); 
void setzd(double val); 
void setk(double val); 
void setzk(double val); 
void setl(double val); 
void setzl(double val); 
void setf(double val); 
void setq(double val); 
void setkci(double val); 
void setkic(double val); 
void setkoi(double val); 
void setkio(double val); 
void setkii2(double val); 
void setki2i(double val); 
void set_log_SGate(bool bLogSGate); 
void set_log_C0Gate(bool bLogC0Gate); 
void set_log_C1Gate(bool bLogC1Gate); 
void set_log_C2Gate(bool bLogC2Gate); 
void set_log_C3Gate(bool bLogC3Gate); 
void set_log_C4Gate(bool bLogC4Gate); 
void set_log_C5Gate(bool bLogC5Gate); 
void set_log_I0Gate(bool bLogI0Gate); 
void set_log_I1Gate(bool bLogI1Gate); 
void set_log_I2Gate(bool bLogI2Gate); 
void set_log_I3Gate(bool bLogI3Gate); 
void set_log_I4Gate(bool bLogI4Gate); 
void set_log_I5Gate(bool bLogI5Gate); 
void set_log_OGate(bool bLogOGate); 
void set_log_I6Gate(bool bLogI6Gate); 
void set_log_I7Gate(bool bLogI7Gate); 

 
protected: 
private: 
 
virtual void specify_write_function_indices(); 
ADouble SGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaSGate; 
ADouble C0Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC0Gate; 
ADouble C1Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC1Gate; 
ADouble C2Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC2Gate; 
ADouble C3Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC3Gate; 
ADouble C4Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC4Gate; 
ADouble C5Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaC5Gate; 
ADouble I0Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI0Gate; 
ADouble I1Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI1Gate; 
ADouble I2Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI2Gate; 
ADouble I3Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI3Gate; 
ADouble I4Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI4Gate; 
ADouble I5Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI5Gate; 
ADouble OGate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaOGate; 
ADouble I6Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI6Gate; 
ADouble I7Gate; 
Grid::AttachmentAccessor<Vertex, ADouble> aaI7Gate; 
number       	gmax ; 
number       	ek; 
number 		a ; 
number 		za ; 
number       	b ; 
number       	zb ; 
number 		c ; 
number 		zc ; 
number 		d ; 
number 		zd ; 
number 		k ; 
number 		zk ; 
number 		l ; 
number 		zl ; 
number 		f ; 
number 		q ; 
number 		kci ; 
number 		kic ; 
number 		koi ; 
number 		kio ; 
number 		kii2 ; 
number 		ki2i ; 
bool m_log_SGate; 
bool m_log_C0Gate; 
bool m_log_C1Gate; 
bool m_log_C2Gate; 
bool m_log_C3Gate; 
bool m_log_C4Gate; 
bool m_log_C5Gate; 
bool m_log_I0Gate; 
bool m_log_I1Gate; 
bool m_log_I2Gate; 
bool m_log_I3Gate; 
bool m_log_I4Gate; 
bool m_log_I5Gate; 
bool m_log_OGate; 
bool m_log_I6Gate; 
bool m_log_I7Gate; 
// Standard-NModl-File-Params 
number F, R, K, celsius; 
}; 
 
} // namespace cable
} // namespace ug


#endif // Kv4_csiosi_converted_standard_UG_H_
