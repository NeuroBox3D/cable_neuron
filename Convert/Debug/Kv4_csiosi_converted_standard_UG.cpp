#include "Kv4_csiosi_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getgmax() 
{ 
return gmax; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getek() 
{ 
return ek; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::geta() 
{ 
return a; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getza() 
{ 
return za; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getb() 
{ 
return b; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getzb() 
{ 
return zb; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getc() 
{ 
return c; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getzc() 
{ 
return zc; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getd() 
{ 
return d; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getzd() 
{ 
return zd; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getk() 
{ 
return k; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getzk() 
{ 
return zk; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getl() 
{ 
return l; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getzl() 
{ 
return zl; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getf() 
{ 
return f; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getq() 
{ 
return q; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getkci() 
{ 
return kci; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getkic() 
{ 
return kic; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getkoi() 
{ 
return koi; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getkio() 
{ 
return kio; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getkii2() 
{ 
return kii2; 
} 
template<typename TDomain> 
double Kv4_csiosi_converted_standard_UG<TDomain>::getki2i() 
{ 
return ki2i; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setgmax(double val) 
{ 
gmax = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setek(double val) 
{ 
ek = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::seta(double val) 
{ 
a = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setza(double val) 
{ 
za = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setb(double val) 
{ 
b = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setzb(double val) 
{ 
zb = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setc(double val) 
{ 
c = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setzc(double val) 
{ 
zc = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setd(double val) 
{ 
d = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setzd(double val) 
{ 
zd = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setk(double val) 
{ 
k = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setzk(double val) 
{ 
zk = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setl(double val) 
{ 
l = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setzl(double val) 
{ 
zl = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setf(double val) 
{ 
f = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setq(double val) 
{ 
q = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setkci(double val) 
{ 
kci = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setkic(double val) 
{ 
kic = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setkoi(double val) 
{ 
koi = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setkio(double val) 
{ 
kio = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setkii2(double val) 
{ 
kii2 = val; 
} 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::setki2i(double val) 
{ 
ki2i = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->C0Gate)) 
UG_THROW("Attachment necessary (C0Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C0Gate); 
this->aaC0Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C0Gate); 
 
if (spGrid->has_vertex_attachment(this->C1Gate)) 
UG_THROW("Attachment necessary (C1Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C1Gate); 
this->aaC1Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C1Gate); 
 
if (spGrid->has_vertex_attachment(this->C2Gate)) 
UG_THROW("Attachment necessary (C2Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C2Gate); 
this->aaC2Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C2Gate); 
 
if (spGrid->has_vertex_attachment(this->C3Gate)) 
UG_THROW("Attachment necessary (C3Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C3Gate); 
this->aaC3Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C3Gate); 
 
if (spGrid->has_vertex_attachment(this->C4Gate)) 
UG_THROW("Attachment necessary (C4Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C4Gate); 
this->aaC4Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C4Gate); 
 
if (spGrid->has_vertex_attachment(this->C5Gate)) 
UG_THROW("Attachment necessary (C5Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->C5Gate); 
this->aaC5Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->C5Gate); 
 
if (spGrid->has_vertex_attachment(this->I0Gate)) 
UG_THROW("Attachment necessary (I0Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I0Gate); 
this->aaI0Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I0Gate); 
 
if (spGrid->has_vertex_attachment(this->I1Gate)) 
UG_THROW("Attachment necessary (I1Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I1Gate); 
this->aaI1Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I1Gate); 
 
if (spGrid->has_vertex_attachment(this->I2Gate)) 
UG_THROW("Attachment necessary (I2Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I2Gate); 
this->aaI2Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I2Gate); 
 
if (spGrid->has_vertex_attachment(this->I3Gate)) 
UG_THROW("Attachment necessary (I3Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I3Gate); 
this->aaI3Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I3Gate); 
 
if (spGrid->has_vertex_attachment(this->I4Gate)) 
UG_THROW("Attachment necessary (I4Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I4Gate); 
this->aaI4Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I4Gate); 
 
if (spGrid->has_vertex_attachment(this->I5Gate)) 
UG_THROW("Attachment necessary (I5Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I5Gate); 
this->aaI5Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I5Gate); 
 
if (spGrid->has_vertex_attachment(this->OGate)) 
UG_THROW("Attachment necessary (OGate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OGate); 
this->aaOGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OGate); 
 
if (spGrid->has_vertex_attachment(this->I6Gate)) 
UG_THROW("Attachment necessary (I6Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I6Gate); 
this->aaI6Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I6Gate); 
 
if (spGrid->has_vertex_attachment(this->I7Gate)) 
UG_THROW("Attachment necessary (I7Gate) for Kv4_csiosi_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->I7Gate); 
this->aaI7Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->I7Gate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> Kv4_csiosi_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
{ 
	 //var for output 
	 std::vector<number> GatingAccesors; 
 
	 typedef ug::MathVector<TDomain::dim> position_type; 
 
	 position_type coord; 
 
	 if (coord.size()==1) 
	 	 coord[0]=x; 
	 if (coord.size()==2) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 } 
	 if (coord.size()==3) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 	 coord[2] = z;
	 } 
	 //accesors 
	 typedef Attachment<position_type> position_attachment_type; 
	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accesor_type; 
 
	 // Definitions for Iteration over all Elements 
	 typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
	 SubsetGroup ssGrp; 
	 try { ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_C0Gate || m_log_C1Gate || m_log_C2Gate || m_log_C3Gate || m_log_C4Gate || m_log_C5Gate || m_log_I0Gate || m_log_I1Gate || m_log_I2Gate || m_log_I3Gate || m_log_I4Gate || m_log_I5Gate || m_log_OGate || m_log_I6Gate || m_log_I7Gate )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pVMDisc->approx_space()->domain()->position_accessor(); 
	 	 	 if (si==0) 
	 	 	 { 
	 	 	 	 bestVrt = *iterBegin; 
	 	 	 	 bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]); 
	 	 	 } 
	 	 	 iter = iterBegin; 
	 	 	 iter++; 
	 	 	 while(iter != iterEnd) 
	 	 	 { 
	 	 	 	 distSq = VecDistanceSq(coord, aaPos[*iter]); 
	 	 	 	 { 
	 	 	 	 	 bestDistSq = distSq; 
	 	 	 	 	 bestVrt = *iter; 
	 	 	 	 } 
	 	 	 	 ++iter; 
	 	 	 } 
	 	 } 
	 	 if (m_log_C0Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC0Gate[bestVrt]); 
	 	 if (m_log_C1Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC1Gate[bestVrt]); 
	 	 if (m_log_C2Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC2Gate[bestVrt]); 
	 	 if (m_log_C3Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC3Gate[bestVrt]); 
	 	 if (m_log_C4Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC4Gate[bestVrt]); 
	 	 if (m_log_C5Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaC5Gate[bestVrt]); 
	 	 if (m_log_I0Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI0Gate[bestVrt]); 
	 	 if (m_log_I1Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI1Gate[bestVrt]); 
	 	 if (m_log_I2Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI2Gate[bestVrt]); 
	 	 if (m_log_I3Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI3Gate[bestVrt]); 
	 	 if (m_log_I4Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI4Gate[bestVrt]); 
	 	 if (m_log_I5Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI5Gate[bestVrt]); 
	 	 if (m_log_OGate == true) 
	 	 	 GatingAccesors.push_back(this->aaOGate[bestVrt]); 
	 	 if (m_log_I6Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI6Gate[bestVrt]); 
	 	 if (m_log_I7Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaI7Gate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C0Gate(bool bLogC0Gate) { m_log_C0Gate = bLogC0Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C1Gate(bool bLogC1Gate) { m_log_C1Gate = bLogC1Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C2Gate(bool bLogC2Gate) { m_log_C2Gate = bLogC2Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C3Gate(bool bLogC3Gate) { m_log_C3Gate = bLogC3Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C4Gate(bool bLogC4Gate) { m_log_C4Gate = bLogC4Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_C5Gate(bool bLogC5Gate) { m_log_C5Gate = bLogC5Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I0Gate(bool bLogI0Gate) { m_log_I0Gate = bLogI0Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I1Gate(bool bLogI1Gate) { m_log_I1Gate = bLogI1Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I2Gate(bool bLogI2Gate) { m_log_I2Gate = bLogI2Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I3Gate(bool bLogI3Gate) { m_log_I3Gate = bLogI3Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I4Gate(bool bLogI4Gate) { m_log_I4Gate = bLogI4Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I5Gate(bool bLogI5Gate) { m_log_I5Gate = bLogI5Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_OGate(bool bLogOGate) { m_log_OGate = bLogOGate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I6Gate(bool bLogI6Gate) { m_log_I6Gate = bLogI6Gate; }
template<typename TDomain> void Kv4_csiosi_converted_standard_UG<TDomain>::set_log_I7Gate(bool bLogI7Gate) { m_log_I7Gate = bLogI7Gate; }
 // Init Method for using gatings 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[VMDisc<TDomain>::_v_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 

 
}  
 
 
 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[VMDisc<TDomain>::_v_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 

 
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
double I6 = aaI6Gate[vrt]; 
double I7 = aaI7Gate[vrt]; 

 
 

 
 
double       kC01f = 4*a*exp(za*v*F/(R*(273.16+celsius)))		;//closed to open pathway transitions; 
double       kC01b = b*exp(zb*v*F/(R*(273.16+celsius))); 
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
double       kOIf = koi						;//open to inactivated transitions		; 
double       kOIb = kio; 
double       kII2f = kii2					; 
double       kII2b = ki2i; 
 
 
 
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
O+=(-O*kOIf+I6*kOIb)*dt; 
I6+=(O*kOIf+-I6*kOIb)*dt; 
I6+=(-I6*kII2f+I7*kII2b)*dt; 
I7+=(I6*kII2f+-I7*kII2b)*dt; 
 
 
 
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
aaI6Gate[vrt] = I6; 
aaI7Gate[vrt] = I7; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
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
number I6 = aaI6Gate[ver]; 
number I7 = aaI7Gate[ver]; 
number k = vrt_values[m_pVMDisc->_k_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 
number g = gmax * O; 

 
 
outCurrentValues.push_back( g * (v - ek)); 
} 
 
 
template<typename TDomain> 
void Kv4_csiosi_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class Kv4_csiosi_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class Kv4_csiosi_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class Kv4_csiosi_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


