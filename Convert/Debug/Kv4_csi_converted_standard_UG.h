#ifndef Kv4_csi_converted_standard_UG_H_
#define Kv4_csi_converted_standard_UG_H_
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
class Kv4_csi_converted_standard_UG
    : public IChannel<TDomain> 
{ 
    public: 
 using IChannel <TDomain>::m_pVMDisc; 
 

 
 



/// @copydoc IChannel<TDomain>::IChannel(cont char*) 
Kv4_csi_converted_standard_UG(const char* functions, const char* subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
      	gmax ( 0.00     	*0.01), 
	ek(-90                  *1), 
	F ( 9.6485e4 			*1), 
      	R ( 8.3145e3 			*1), 
	a ( 7			*1), 
	za ( 0.315646648				*1), 
      	b ( .090		*1), 
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
	kic ( 0.0017935468 	*1) {} 
UG_CATCH_THROW("Error in Kv4_csi_converted_standard_UG initializer list. ") 
 
 
/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&) 
Kv4_csi_converted_standard_UG(const std::vector<std::string>& functions, const std::vector<std::string>& subsets) 
try : IChannel<TDomain>(functions, subsets), 
m_R(8.314), m_T(293.0), m_F(96485.0), 
      	gmax ( 0.00     	*0.01), 
	ek(-90                  *1), 
	F ( 9.6485e4 			*1), 
      	R ( 8.3145e3 			*1), 
	a ( 7			*1), 
	za ( 0.315646648				*1), 
      	b ( .090		*1), 
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
	kic ( 0.0017935468 	*1) {} 
UG_CATCH_THROW("Error in Kv4_csi_converted_standard_UG initializer list. ") 
/// destructor 
 
virtual ~Kv4_csi_converted_standard_UG() {}; 
/// create attachments and accessors 
void init_attachments(); 
// inherited from IChannel 
 
virtual void init(const LocalVector& u, Edge* e); 
virtual void update_gating(number newTime, const LocalVector& u, Edge* e); 
virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues); 
virtual void vm_disc_available(); 

 
double getgmax(); 
double getek(); 
double getF(); 
double getR(); 
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
void setgmax(double val); 
void setek(double val); 
void setF(double val); 
void setR(double val); 
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

 
protected: 
private: 
 
number m_R, m_T, m_F; 
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
number       	gmax ; 
number 	ek; 
number 	F ; 
number       	R ; 
number 	a ; 
number 	za ; 
number       	b ; 
number       	zb ; 
number 	c ; 
number 	zc ; 
number 	d ; 
number 	zd ; 
number 	k ; 
number 	zk ; 
number 	l ; 
number 	zl ; 
number 	f ; 
number 	q ; 
number 	kci ; 
number 	kic ; 
}; 
 
} // namespace cable
} // namespace ug
 
 
#endif // Kv4_csi_converted_standard_UG_H_
