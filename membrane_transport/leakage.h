/*
 * leakage.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__


#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable {


template <typename TDomain>
class ChannelLeak
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		ChannelLeak(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		ChannelLeak(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~ChannelLeak() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set leakage conductance
		void set_cond(number g);
		/// set leakage conductance on specific subsets
		/// \{
		void set_cond(number g, const char* subsets);
		void set_cond(number g, const std::vector<std::string>& subsets);
		/// \}

		/// set leakage equilibrium potential
		void set_rev_pot(number e);
		/// set leakage equilibrium potential on specific subsets
		void set_rev_pot(number e, const char* subsets);
		void set_rev_pot(number e, const std::vector<std::string>& subsets);


		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z);
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		struct Params
		{
			Params() : g(1.0e-6), E(-65.0) {};
			number g; ///< membrane conductance [C / (m^2 * mV * ms)]
			number E; ///< reversal potential [mV]
		};

		std::map<std::string, Params> m_mSubsetParams2Save;
		std::map<int, Params> m_mSubsetParams;
};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__LEAKAGE_H__
