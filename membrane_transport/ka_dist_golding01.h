/*
 * kadist_g01.h
 *
 *  Created on: 2019-05-17
 *      Author: mbreit
 */


#ifndef __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KADIST_G01_H__
#define __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KADIST_G01_H__

#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


// forward declaration 
template <typename TDomain> 
class CableEquation; 


/**
 * @brief K channels from Golding 2001
 */
template <typename TDomain> 
class KaDist_Golding01
: public ICableMembraneTransport<TDomain>
{ 
	public:
		using ICableMembraneTransport <TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*, const char*)
		KaDist_Golding01(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&, const std::vector<std::string>&)
		KaDist_Golding01
		(
			const std::vector<std::string>& functions,
			const std::vector<std::string>& subsets
		);

		/// destructor
		virtual ~KaDist_Golding01();
		// inherited from ICableMembraneTransport

		// setters
		void set_ek(number val);
		void set_gkabar(number val);
		void set_vhalfn(number val);
		void set_vhalfl(number val);
		void set_a0n(number val);
		void set_zetan(number val);
		void set_zetal(number val);
		void set_gmn(number val);
		void set_gml(number val);
		void set_lmin(number val);
		void set_nmin(number val);
		void set_pw(number val);
		void set_tq(number val);
		void set_qq(number val);
		void set_q10(number val);

		void set_logLGate(bool bLoglGate);
		void set_logNGate(bool bLognGate);


		// getters
		number ek() const;
		number gkabar() const;
		number vhalfn() const;
		number vhalfl() const;
		number a0n() const;
		number zetan() const;
		number zetal() const;
		number gmn() const;
		number gml() const;
		number lmin() const;
		number nmin() const;
		number pw() const;
		number tq() const;
		number qq() const;
		number q10() const;

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
		number alpl(number v);
		number betl(number v);

		/// create attachments and accessors
		void init_attachments();


	private:
		ANumber m_aNGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaNGate;
		ANumber m_aLGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaLGate;

		number m_ek;  // V
		number m_gkabar;  // S/m^2

		number m_vhalfn;  // V
		number m_vhalfl;  // V
		//number a0l;
		number m_a0n;  // s^-1
		number m_zetan;  // 1
		number m_zetal;  // 1
		number m_gmn;  // 1
		number m_gml;  // 1
		number m_lmin;  // s
		number m_nmin;  // s
		number m_pw;  // 1
		number m_tq;  // V
		number m_qq;  // V

		number m_q10;  // 1
		//number qtl;

		bool m_bLogNGate;
		bool m_bLogLGate;
}; 

} // namespace cable_neuron
} // namespace ug


#endif // __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__KADIST_G01_H__
