#include "cad_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getdepth() 
{ 
return depth; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::gettaur() 
{ 
return taur; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getcainf() 
{ 
return cainf; 
} 
template<typename TDomain> 
double cad_converted_standard_UG<TDomain>::getcai() 
{ 
return cai; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setdepth(double val) 
{ 
depth = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::settaur(double val) 
{ 
taur = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setcainf(double val) 
{ 
cainf = val; 
} 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::setcai(double val) 
{ 
cai = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->caSGate)) 
UG_THROW("Attachment necessary (caSGate) for cad_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->caSGate); 
this->aacaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->caSGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> cad_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_caSGate == true )
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
	 	 if (m_log_caSGate == true) 
	 	 	 GatingAccesors.push_back(this->aacaSGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void cad_converted_standard_UG<TDomain>::set_log_caSGate(bool bLogcaSGate) { m_log_caSGate = bLogcaSGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
number dt = m_pVMDisc->time(); 
number ica = m_pVMDisc->flux_ca(); 
// make preparing vor getting values of every edge 
number v = vrt_values[VMDisc<TDomain>::_v_]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 

 
aacaSGate[vrt] =  cainf; 
cai =  ca; 
}  
 
 
 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number ica = m_pVMDisc->flux_ca(); 
number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[VMDisc<TDomain>::_v_]; 
number ca = vrt_values[VMDisc<TDomain>::_ca_]; 

 
double caS = aacaSGate[vrt]; 

 
 
double 	drive_channel =  - (10000) * ica / (2 * FARADAY * depth); 
; 
 
if (drive_channel <= 0.)
{ 
 drive_channel = 0. ; 
} 
caS +=  drive_channel + (cainf-caS)/taur*dt; 
; 
 
	cai = ca; 
; 
 

 
 
aacaSGate[vrt] = caS; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void cad_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number ica = m_pVMDisc->flux_ca(); 
number caS = aacaSGate[ver]; 
number ca = vrt_values[m_pVMDisc->_ca_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 

 
 
outCurrentValues.push_back(   caS); 
outCurrentValues.push_back(   caS ); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class cad_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class cad_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class cad_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


