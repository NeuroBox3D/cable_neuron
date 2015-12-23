#include "caL3d_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::ghk(double v, double  ci, double  co) 
{ 

double 	z = (0.001)*2*F*v/(R*(celsius+273.15)); 
	return  (.001)*2*F*(ci*efun(-z) - co*efun(z)); 
}
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::efun(double z) 
{ 
double efun; 
	if (fabs(z) < 1e-4) {
		return  1 - z/2; 
	}else{
		return  z/(exp(z) - 1); 
	}
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getp() 
{ 
return p; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getth() 
{ 
return th; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getq() 
{ 
return q; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getRa() 
{ 
return Ra; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::gettemp() 
{ 
return temp; 
} 
template<typename TDomain> 
double caL3d_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setp(double val) 
{ 
p = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setth(double val) 
{ 
th = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setq(double val) 
{ 
q = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setRa(double val) 
{ 
Ra = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::settemp(double val) 
{ 
temp = val; 
} 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->CGate)) 
UG_THROW("Attachment necessary (CGate) for caL3d_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->CGate); 
this->aaCGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->CGate); 
 
if (spGrid->has_vertex_attachment(this->OGate)) 
UG_THROW("Attachment necessary (OGate) for caL3d_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OGate); 
this->aaOGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> caL3d_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_CGate || m_log_OGate )
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
	 	 if (m_log_CGate == true) 
	 	 	 GatingAccesors.push_back(this->aaCGate[bestVrt]); 
	 	 if (m_log_OGate == true) 
	 	 	 GatingAccesors.push_back(this->aaOGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void caL3d_converted_standard_UG<TDomain>::set_log_CGate(bool bLogCGate) { m_log_CGate = bLogCGate; }
template<typename TDomain> void caL3d_converted_standard_UG<TDomain>::set_log_OGate(bool bLogOGate) { m_log_OGate = bLogOGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 

 
}  
 
 
 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
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
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 

 
double C = aaCGate[vrt]; 
double O = aaOGate[vrt]; 

 
 

 
 
tadj= pow(q10 , ((celsius-temp)/10)); 
number 	a = Ra / (1 + exp(-(v-th)/q)) * tadj; 
number 	b = Rb / (1 + exp((v-th)/q)) * tadj; 
 
 
 
C+=(-C*a+O*b)*dt; 
O+=(C*a+-O*b)*dt; 
 
 
 
aaCGate[vrt] = C; 
aaOGate[vrt] = O; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number C = aaCGate[ver]; 
number O = aaOGate[ver]; 
number ca = vrt_values[m_pVMDisc->_ca_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
number cai = ca;

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 
number cao = m_pVMDisc->ca_out(); 

 
 

 
 
number rates(v); 
outCurrentValues.push_back( O * p * ghk(v,cai,cao)); 
} 
 
 
template<typename TDomain> 
void caL3d_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class caL3d_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class caL3d_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class caL3d_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


