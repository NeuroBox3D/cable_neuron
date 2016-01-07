#include "../converted/kir_converted_standard_UG.h"
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::taumkir(double v)
{  
return v; 
} 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::ce_obj_available()  
{  
	init_attachments();  
 	F = m_pCE->F; 
 R = m_pCE->R; 
 K = m_pCE->temperature(); 
 celsius = m_pCE->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getgkbar() 
{ 
return gkbar; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmvhalf() 
{ 
return mvhalf; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmslope() 
{ 
return mslope; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getmshift() 
{ 
return mshift; 
} 
template<typename TDomain> 
double kir_converted_standard_UG<TDomain>::getqfact() 
{ 
return qfact; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setgkbar(double val) 
{ 
gkbar = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmvhalf(double val) 
{ 
mvhalf = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmslope(double val) 
{ 
mslope = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setmshift(double val) 
{ 
mshift = val; 
} 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::setqfact(double val) 
{ 
qfact = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for kir_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> kir_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 try { ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_mGate )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pCE->approx_space()->domain()->position_accessor(); 
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
	 	 if (m_log_mGate == true) 
	 	 	 GatingAccesors.push_back(this->aamGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void kir_converted_standard_UG<TDomain>::set_log_mGate(bool bLogmGate) { m_log_mGate = bLogmGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
number dt = m_pCE->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
double 			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) ); 
aamGate[vrt] = minf; 
}  
 
 
 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
 number FARADAY = m_pCE->F; 
 number dt = newTime - m_pCE->time(); 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 

 
double m = aamGate[vrt]; 

 
 
double 			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) ); 
        m  +=  (minf - m) / ( taumkir(v)/qfact )*dt; 
; 
 

 
 
aamGate[vrt] = m; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number m = aamGate[ver]; 
number k = vrt_values[m_pCE->_k_]; 
number v =  vrt_values[m_pCE->_v_]; 
 
 
number t = m_pCE->time(); 
 
 

 
 
const number helpV = 1e3*(m_pCE->R*m_pCE->temperature())/m_pCE->F; 
number ek; 
if (m_pCE->rev_pot_k() == 0) 
{ 
	  ek = helpV*(log(m_pCE->k_out()/k)); 
} 
else 
{ 
	  ek = m_pCE->rev_pot_k(); 
} 
 
 
number gk = gkbar * m; 

 
 
outCurrentValues.push_back( gk * ( v - ek )); 
} 
 
 
template<typename TDomain> 
void kir_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class kir_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class kir_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class kir_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


