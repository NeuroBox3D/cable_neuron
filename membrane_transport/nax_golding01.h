/*
 * nax_golding01.h
 *
 *  Created on: 2019-06-12
 *      Author: mbreit
 */


#ifndef UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H__
#define UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H__

#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 


/**
 * @brief Na channel for axons
 * without slow inactivation
 * taken from Golding 2001
 */
template <typename TDomain>
class Nax_Golding01
: public ICableMembraneTransport<TDomain>
{ 
	public:
		using ICableMembraneTransport <TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*, const char*)
		Nax_Golding01(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&, const std::vector<std::string>&)
		Nax_Golding01
		(
			const std::vector<std::string>& functions,
			const std::vector<std::string>& subsets
		);

		/// destructor
		virtual ~Nax_Golding01();

		/// name
		virtual std::string name();


		// setters
		void set_gbar(number val);
		void set_tha(number val);
		void set_qa(number val);
		void set_Ra(number val);
		void set_Rb(number val);
		void set_thi1(number val);
		void set_thi2(number val);
		void set_qd(number val);
		void set_qg(number val);
		void set_mmin(number val);
		void set_hmin(number val);
		void set_q10(number val);
		void set_Rg(number val);
		void set_Rd(number val);
		void set_thinf(number val);
		void set_qinf(number val4);

#ifdef UG_FOR_LUA
		// function setter for conductance
		void set_gbar(SmartPtr<LuaUserData<number, TDomain::dim> > fct);
		void set_gbar(const char* fct);
#endif

		void set_logMGate(bool b);
		void set_logHGate(bool b);


		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newtime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;

	private:
		virtual void specify_write_function_indices();

	protected:
		number alpn(number v);
		number betn(number v);

		/// create attachments and accessors
		void init_attachments();


	private:
		ANumber m_aMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		ANumber m_aHGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;

#ifdef UG_FOR_LUA
		bool m_bConductanceDependsOnCoordinates;
		SmartPtr<LuaUserData<number, TDomain::dim> > m_spCondFct;
		ANumber m_aGBar;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaGBar;
#endif
		number m_gbar;  // S/m^2

		number m_tha;  // V
		number m_qa;  // V
		number m_Ra;  // s^-1
		number m_Rb;  // s^-1

		number m_thi1;  // V
		number m_thi2;  // V
		number m_qd;  // V
		number m_qg;  // V

		number m_mmin;  // s
		number m_hmin;  // s
		number m_q10;  // 1
		number m_Rg;  // s^-1
		number m_Rd;  // s^-1

		number m_thinf;  // V
		number m_qinf;  // V

		bool m_bLogMGate;
		bool m_bLogHGate;
}; 

} // namespace cable_neuron
} // namespace ug


#endif // UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__NAX_GOLDING01_H__
