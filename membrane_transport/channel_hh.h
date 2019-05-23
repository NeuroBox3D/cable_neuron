/*
 * channel_hh.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CHANNEL_HH_H__
#define __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CHANNEL_HH_H__


#include "cable_membrane_transport_interface.h"

namespace ug {
namespace cable_neuron {


template <typename TDomain>
class ChannelHH
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		ChannelHH(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		ChannelHH(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);;

		/// destructor
		virtual ~ChannelHH() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set conductances
		void set_conductances(number gK, number gNa);
		/// set conductances on specific subsets
		/// \{
		void set_conductances(number gK, number gNa, const char* subsets);
		void set_conductances(number gK, number gNa, const std::vector<std::string>& subsets);
		/// \}

		/// whether or not temperature dependency is to be enabled
		void enable_temperature_dependency(bool enable);

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	private:
		number vtrap(number x, number y);

	public:
		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		struct Params
		{
			Params() : gK(3.6e2), gNa(1.2e3) {};
			number gK;	///< potassium conductance [S/m^2)]
			number gNa;	///< sodium conductance [S/m^2]
		};

		std::map<std::string, Params> m_mSubsetParams2Save;
		std::map<int, Params> m_mSubsetParams;

		// one attachment per state variable
		ANumber m_MGate;
		ANumber m_HGate;
		ANumber m_NGate;

		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaNGate;

		// attachment log_files
		bool m_log_nGate, m_log_hGate, m_log_mGate;

		// temperature dependency flag
		bool m_bTempDep;
};



template <typename TDomain>
class ChannelHHNernst
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		ChannelHHNernst(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		ChannelHHNernst(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~ChannelHHNernst() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set conductivities
		void set_conductances(number gK, number gNa);
		/// set conductances on specific subsets
		/// \{
		void set_conductances(number gK, number gNa, const char* subsets);
		void set_conductances(number gK, number gNa, const std::vector<std::string>& subsets);
		/// \}

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	private:
		number vtrap(number x, number y);

	public:
		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z) const;
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		struct Params
		{
			Params() : gK(3.6e2), gNa(1.2e3) {};
			number gK;	///< potassium conductance [S/m^2)]
			number gNa;	///< sodium conductance [S/m^2]
		};

		std::map<std::string, Params> m_mSubsetParams2Save;
		std::map<int, Params> m_mSubsetParams;

		// one attachment per state variable
		ANumber m_MGate;
		ANumber m_HGate;
		ANumber m_NGate;

		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaNGate;

		// attachment log_files
		bool m_log_nGate, m_log_hGate, m_log_mGate;
};

} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__CHANNEL_HH_H__
