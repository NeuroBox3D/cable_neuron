#include "Kv4_csi_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable {
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getgmax() 
{ 
return gmax; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getek() 
{ 
return ek; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getF() 
{ 
return F; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getR() 
{ 
return R; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::geta() 
{ 
return a; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getza() 
{ 
return za; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getb() 
{ 
return b; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getzb() 
{ 
return zb; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getc() 
{ 
return c; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getzc() 
{ 
return zc; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getd() 
{ 
return d; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getzd() 
{ 
return zd; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getk() 
{ 
return k; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getzk() 
{ 
return zk; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getl() 
{ 
return l; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getzl() 
{ 
return zl; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getf() 
{ 
return f; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getq() 
{ 
return q; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getkci() 
{ 
return kci; 
} 
template<typename TDomain> 
double Kv4_csi_converted_standard_UG<TDomain>::getkic() 
{ 
return kic; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setgmax(double val) 
{ 
gmax = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setek(double val) 
{ 
ek = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setF(double val) 
{ 
F = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setR(double val) 
{ 
R = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::seta(double val) 
{ 
a = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setza(double val) 
{ 
za = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setb(double val) 
{ 
b = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setzb(double val) 
{ 
zb = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setc(double val) 
{ 
c = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setzc(double val) 
{ 
zc = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setd(double val) 
{ 
d = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setzd(double val) 
{ 
zd = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setk(double val) 
{ 
k = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setzk(double val) 
{ 
zk = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setl(double val) 
{ 
l = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setzl(double val) 
{ 
zl = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setf(double val) 
{ 
f = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setq(double val) 
{ 
q = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setkci(double val) 
{ 
kci = val; 
} 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::setkic(double val) 
{ 
kic = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->C0Gate)) 
UG_THROW("Attachment necessary (C0Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C0Gate); 
this->aaC0Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C0Gate); 
 
if (spGrid->has_vertex_attachment(this->C1Gate)) 
UG_THROW("Attachment necessary (C1Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C1Gate); 
this->aaC1Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C1Gate); 
 
if (spGrid->has_vertex_attachment(this->C2Gate)) 
UG_THROW("Attachment necessary (C2Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C2Gate); 
this->aaC2Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C2Gate); 
 
if (spGrid->has_vertex_attachment(this->C3Gate)) 
UG_THROW("Attachment necessary (C3Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C3Gate); 
this->aaC3Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C3Gate); 
 
if (spGrid->has_vertex_attachment(this->C4Gate)) 
UG_THROW("Attachment necessary (C4Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C4Gate); 
this->aaC4Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C4Gate); 
 
if (spGrid->has_vertex_attachment(this->C5Gate)) 
UG_THROW("Attachment necessary (C5Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C5Gate); 
this->aaC5Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C5Gate); 
 
if (spGrid->has_vertex_attachment(this->I0Gate)) 
UG_THROW("Attachment necessary (I0Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I0Gate); 
this->aaI0Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I0Gate); 
 
if (spGrid->has_vertex_attachment(this->I1Gate)) 
UG_THROW("Attachment necessary (I1Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I1Gate); 
this->aaI1Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I1Gate); 
 
if (spGrid->has_vertex_attachment(this->I2Gate)) 
UG_THROW("Attachment necessary (I2Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I2Gate); 
this->aaI2Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I2Gate); 
 
if (spGrid->has_vertex_attachment(this->I3Gate)) 
UG_THROW("Attachment necessary (I3Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I3Gate); 
this->aaI3Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I3Gate); 
 
if (spGrid->has_vertex_attachment(this->I4Gate)) 
UG_THROW("Attachment necessary (I4Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I4Gate); 
this->aaI4Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I4Gate); 
 
if (spGrid->has_vertex_attachment(this->I5Gate)) 
UG_THROW("Attachment necessary (I5Gate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I5Gate); 
this->aaI5Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I5Gate); 
 
if (spGrid->has_vertex_attachment(this->OGate)) 
UG_THROW("Attachment necessary (OGate) for Kv4_csi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OGate); 
this->aaOGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius and time
number celsius = m_pVMDisc->celsius; 
number dt = m_pVMDisc->time(); 
number ik = m_pVMDisc->get_flux_k(); 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number v = u(m_pVMDisc->_v_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
}  
}  
 
 
 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pVMDisc->celsius; 
 number FARADAY = m_F; 
 number ik = m_pVMDisc->get_flux_k(); 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number dt = newTime - m_pVMDisc->m_aaTime[vrt]; 
number v = u(m_pVMDisc->_v_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
double C0 = aaC0Gate[vrt]; 
double C1 = aaC1Gate[vrt]; 
double C2 = aaC2Gate[vrt]; 
double C3 = aaC3Gate[vrt]; 
double C4 = aaC4Gate[vrt]; 
double C5 = aaC5Gate[vrt]; 
double I0 = aaI0Gate[vrt]; 
double I1 = aaI1Gate[vrt]; 
double I2 = aaI2Gate[vrt]; 
double I3 = aaI3Gate[vrt]; 
double I4 = aaI4Gate[vrt]; 
double I5 = aaI5Gate[vrt]; 
double O = aaOGate[vrt]; 

 
 

 
 
double       kC01f = 4*a*exp(za*v*F/(R*(273.16+celsius)))		;//closed to open pathway transitions; 
double       kC01b = b*exp(zb*v*F/(R*(273.16+celsius)))		;//273.16 K = 0.01degCelsius; 
double       kC12f = 3*a*exp(za*v*F/(R*(273.16+celsius))); 
double       kC12b = 2*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kC23f = 2*a*exp(za*v*F/(R*(273.16+celsius))); 
double       kC23b = 3*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kC34f = a*exp(za*v*F/(R*(273.16+celsius))); 
double       kC34b = 4*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kC45f = c*exp(zc*v*F/(R*(273.16+celsius))); 
double       kC45b = d*exp(zd*v*F/(R*(273.16+celsius))); 
double       kCOf = k*exp(zk*v*F/(R*(273.16+celsius))); 
double       kCOb = l*exp(zl*v*F/(R*(273.16+celsius))); 
double kCI0f=kci*( pow(f , 4));//closedtoinactivatedtransitions; 
double kCI0b=kic/( pow(f , 4)); 
double kCI1f=kci*( pow(f , 3)); 
double kCI1b=kic/( pow(f , 3)); 
double kCI2f=kci*( pow(f , 2)); 
double kCI2b=kic/( pow(f , 2)); 
double       kCI3f = kci*(f); 
double       kCI3b = kic/(f); 
double       kCI4f = kci; 
double       kCI4b = kic; 
double       kCI5f = kci*q; 
double       kCI5b = kic/q; 
double       kI01f = 4*(1/f)*a*exp(za*v*F/(R*(273.16+celsius)))	;//closed to inactivated transitions; 
double       kI01b = d*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kI12f = 3*(1/f)*a*exp(za*v*F/(R*(273.16+celsius))); 
double       kI12b = 2*f*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kI23f = 2*(1/f)*a*exp(za*v*F/(R*(273.16+celsius))); 
double       kI23b = 3*f*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kI34f = (1/f)*a*exp(za*v*F/(R*(273.16+celsius))); 
double       kI34b = 4*f*b*exp(zb*v*F/(R*(273.16+celsius))); 
double       kI45f = q*c*exp(zc*v*F/(R*(273.16+celsius))); 
double       kI45b = (1/q)*d*exp(zd*v*F/(R*(273.16+celsius))); 
 
 
 
C0+=(-C0*kC01f+C1*kC01b)*dt; 
C1+=(C0*kC01f+-C1*kC01b)*dt; 
C1+=(-C1*kC12f+C2*kC12b)*dt; 
C2+=(C1*kC12f+-C2*kC12b)*dt; 
C2+=(-C2*kC23f+C3*kC23b)*dt; 
C3+=(C2*kC23f+-C3*kC23b)*dt; 
C3+=(-C3*kC34f+C4*kC34b)*dt; 
C4+=(C3*kC34f+-C4*kC34b)*dt; 
C4+=(-C4*kC45f+C5*kC45b)*dt; 
C5+=(C4*kC45f+-C5*kC45b)*dt; 
C5+=(-C5*kCOf+O*kCOb)*dt; 
O+=(C5*kCOf+-O*kCOb)*dt; 
C0+=(-C0*kCI0f+I0*kCI0b)*dt; 
I0+=(C0*kCI0f+-I0*kCI0b)*dt; 
C1+=(-C1*kCI1f+I1*kCI1b)*dt; 
I1+=(C1*kCI1f+-I1*kCI1b)*dt; 
C2+=(-C2*kCI2f+I2*kCI2b)*dt; 
I2+=(C2*kCI2f+-I2*kCI2b)*dt; 
C3+=(-C3*kCI3f+I3*kCI3b)*dt; 
I3+=(C3*kCI3f+-I3*kCI3b)*dt; 
C4+=(-C4*kCI4f+I4*kCI4b)*dt; 
I4+=(C4*kCI4f+-I4*kCI4b)*dt; 
C5+=(-C5*kCI5f+I5*kCI5b)*dt; 
I5+=(C5*kCI5f+-I5*kCI5b)*dt; 
I0+=(-I0*kI01f+I1*kI01b)*dt; 
I1+=(I0*kI01f+-I1*kI01b)*dt; 
I1+=(-I1*kI12f+I2*kI12b)*dt; 
I2+=(I1*kI12f+-I2*kI12b)*dt; 
I2+=(-I2*kI23f+I3*kI23b)*dt; 
I3+=(I2*kI23f+-I3*kI23b)*dt; 
I3+=(-I3*kI34f+I4*kI34b)*dt; 
I4+=(I3*kI34f+-I4*kI34b)*dt; 
I4+=(-I4*kI45f+I5*kI45b)*dt; 
I5+=(I4*kI45f+-I5*kI45b)*dt; 
 
 
 
aaC0Gate[vrt] = C0; 
aaC1Gate[vrt] = C1; 
aaC2Gate[vrt] = C2; 
aaC3Gate[vrt] = C3; 
aaC4Gate[vrt] = C4; 
aaC5Gate[vrt] = C5; 
aaI0Gate[vrt] = I0; 
aaI1Gate[vrt] = I1; 
aaI2Gate[vrt] = I2; 
aaI3Gate[vrt] = I3; 
aaI4Gate[vrt] = I4; 
aaI5Gate[vrt] = I5; 
aaOGate[vrt] = O; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void Kv4_csi_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number ik = m_pVMDisc->get_flux_k(); 
number C0 = aaC0Gate[ver]; 
number C1 = aaC1Gate[ver]; 
number C2 = aaC2Gate[ver]; 
number C3 = aaC3Gate[ver]; 
number C4 = aaC4Gate[ver]; 
number C5 = aaC5Gate[ver]; 
number I0 = aaI0Gate[ver]; 
number I1 = aaI1Gate[ver]; 
number I2 = aaI2Gate[ver]; 
number I3 = aaI3Gate[ver]; 
number I4 = aaI4Gate[ver]; 
number I5 = aaI5Gate[ver]; 
number O = aaOGate[ver]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
 
 
number g = gmax * O; 

 
 
number ko = m_pVMDisc->k_out; 

 
 
outCurrentValues.push_back( g * (v - ek)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class Kv4_csi_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class Kv4_csi_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class Kv4_csi_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug
  
  
