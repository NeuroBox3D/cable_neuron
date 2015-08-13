/*
 * channel_interface.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__


#include "VM_Disc.h"


namespace ug {
namespace cable {


// forward declaration
template <typename TDomain>
class VMDisc;


template <typename TDomain>
class IChannel
{
	public:
		// TODO: We do not need functions here!
		///	constructor with comma-separated c-string
		IChannel(const char* functions, const char* subsets);

		/// constructor with vector of string
		IChannel(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		///	destructor
		virtual ~IChannel() {};

		/// name
		virtual std::string name() {return std::string("");};

		/**
		 * @brief Initializes the defined channel type.
		 *
		 * During the initialization, the necessary attachments are attached to the vertices
		 * and their values calculated by the equilibrium state for the initial membrane potential
		 * given in the grid function passed
		 *
		 * @param time			initial time
		 * @param spGridFct		initial solution (containing membrane potential and ion concentrations)
		 */
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// updates the gating parameters
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// provides the ionic current (mol*s^-1) at a given vertex
		virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) = 0;

		/// called when approximation space is available
		void approx_space_available();

		/// called when access to the root VM disc is possible
		virtual void vm_disc_available() {};

		/// getting values of internal channel states
		virtual std::vector<number> state_values(number x, number y, number z) = 0;

		/// adding some Jacobian infos at given vertex
		//virtual void Jacobi_sets(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outJFlux) = 0;

		/**
		 * 	@brief calculate (estimate) linear dependency on potential
		 *
		 *	This method is useful for automatic time step size calculation. It is supposed to calculate
		 *	(or at least estimate) the linear dependency (first term in the taylor expansion)
		 *	of the electric current through the mechanism on the given potential. This is used in
		 *	the VMDisc to get an estimate for the CFL condition that has to be fulfilled by the
		 *	time step size.
		 *
		 * @param vrt			vertex to compute dependency for
		 * @param vrt_values	current solution at this vertex
		 * @return				linear dependency
		 */
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values) {return 0.0;}


		const std::vector<size_t>& fct_indices() const {return m_vWFctInd;}
		const std::vector<std::string>& write_subsets() {return m_vSubset;}

		/**
		 * @brief check definition on a given subset
		 *
		 * @param	si subset index to compare against
		 * @return	true iff channel is defined on subset index si
		 */
		bool is_def_on_subset(int si) const;

		void set_vm_disc(VMDisc<TDomain>*  vmdisc) {m_pVMDisc = vmdisc;}

	private:
		virtual void specify_write_function_indices() {UG_THROW("specify_write_function_indices() not implemented!");}

	protected:
		/// indices in vmdisc for functions whose defect will be written to by this channel
		std::vector<size_t> m_vWFctInd;

		/// vector of subsets this channel is declared on
		std::vector<std::string> m_vSubset;

		/// vector of subset indices this channel is declared on (sorted)
		std::vector<int> m_vSI;

		/// joint VMDisc
		VMDisc<TDomain>* m_pVMDisc;

};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
