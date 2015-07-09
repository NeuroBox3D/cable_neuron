#include "release_exp_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
}  
 
 
 
template<typename TDomain> 
double release_exp_converted_standard_UG<TDomain>::gettau1() 
{ 
return tau1; 
} 
template<typename TDomain> 
double release_exp_converted_standard_UG<TDomain>::gettau2() 
{ 
return tau2; 
} 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::settau1(double val) 
{ 
tau1 = val; 
} 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::settau2(double val) 
{ 
tau2 = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->AGate)) 
UG_THROW("Attachment necessary (AGate) for release_exp_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->AGate); 
this->aaAGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->AGate); 
 
if (spGrid->has_vertex_attachment(this->BGate)) 
UG_THROW("Attachment necessary (BGate) for release_exp_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->BGate); 
this->aaBGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->BGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> release_exp_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_AGate == true || m_log_BGate == true )
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
	 	 if (m_log_AGate == true) 
	 	 	 GatingAccesors.push_back(this->aaAGate[bestVrt]); 
	 	 if (m_log_BGate == true) 
	 	 	 GatingAccesors.push_back(this->aaBGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void release_exp_converted_standard_UG<TDomain>::set_log_AGate(bool bLogAGate) { m_log_AGate = bLogAGate; }
template<typename TDomain> void release_exp_converted_standard_UG<TDomain>::set_log_BGate(bool bLogBGate) { m_log_BGate = bLogBGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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

 
total =  0; 
tau1 =  .9999*tau2; 
aaAGate[vrt] = 0; 
aaBGate[vrt] = 0; 
double tp =  (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1); 
double factor =  -exp(-tp/tau1) + exp(-tp/tau2); 
factor =  1/factor;
}  
 
 
 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
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

 
double A = aaAGate[vrt]; 
double B = aaBGate[vrt]; 

 
 
    B  +=  -B/tau2*dt; 
; 
 

 
 
aaAGate[vrt] = A; 
aaBGate[vrt] = B; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void release_exp_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number A = aaAGate[ver]; 
number B = aaBGate[ver]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 

 
 
number T = B - A; 
outCurrentValues.push_back(0); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class release_exp_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class release_exp_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class release_exp_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


