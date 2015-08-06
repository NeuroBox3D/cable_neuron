/*
 * channel_hh.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_HH_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_HH_H__


#include "channel_interface.h"

namespace ug {
namespace cable {


template <typename TDomain>
class ChannelHH
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		ChannelHH(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		ChannelHH(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);;

		/// destructor
		virtual ~ChannelHH() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set conductivities
		void set_conductances(number gK, number gNa);

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	private:
		number vtrap(number x, number y);

	public:
		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		virtual std::vector<number> state_values(number x, number y, number z);
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		number m_g_K;		// C / (m^2 * mV * ms)
		number m_g_Na;		// C / (m^2 * mV * ms)

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



template <typename TDomain>
class ChannelHHNernst
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		ChannelHHNernst(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		ChannelHHNernst(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~ChannelHHNernst() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// set conductivities
		void set_conductances(number gK, number gNa);

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	private:
		number vtrap(number x, number y);

	public:
		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vmDisc_available();
		virtual std::vector<number> state_values(number x, number y, number z);
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		number m_g_K;
		number m_g_Na;

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

} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_HH_H__
