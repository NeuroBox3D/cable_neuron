/*!
 * \file synapse_info.h
 * \brief synapse info struct
 *
 *  Created on: May 11, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO__
#define __H__UG__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO__

/// includes
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <common/types.h>
#include <common/error.h>	// UG_THROW
#include <common/log.h>		// UG_LOG

#include "../function/types.h"
#include "../function/synapses.h" // for synapse types

/* \defgroup sh_plugin Synapse Handler
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
namespace synapse_handler {


// forward declaration
template <typename TSynapse>
struct synapse_traits;


/*!
 * \brief stores synapse information
 */
struct SynapseInfo {
   /////////////////////////////////////////////
   /// members
   /////////////////////////////////////////////
   number m_locCoords; //!< local coordinates on edge
   unsigned char m_type; //!< type of synapse (\see types.h for types)
   number m_onset; //!< onset of synapse
   number m_tau; //!< decay constant or rise constant
   number m_gMax; //!< g_max - maximal conductance (including all NEURON factors)
   number m_vRev; //!< v_rev - reversal potential (mV)
   number m_param1; //!< generic parameter #1
   number m_param2; //!< generic parameter #2
   number m_param3; //!< generic parameter #3
   unsigned int m_param4; //!< generic parameter #4

  /*!
   * \brief def ctor
   */
   SynapseInfo() :
		   m_locCoords(0), m_type(EMPTY_SYNAPSE), m_onset(std::numeric_limits<number>::quiet_NaN()),
		   m_tau(0), m_gMax(0), m_vRev(0), m_param1(0), m_param2(0), m_param3(0), m_param4(0) { }

   /*!
	* \brief full ctor
	*/
   SynapseInfo(number locCoords, unsigned char type, number onset=std::numeric_limits<number>::quiet_NaN(),
		   number tau=0, number gmax=0, number vrev=0, number param1=0, number param2=0, number param3=0,
		   unsigned int param4=0) :
	   m_locCoords(locCoords), m_type(type), m_onset(onset), m_tau(tau),
	   m_gMax(gmax), m_vRev(vrev), m_param1(param1), m_param2(param2), m_param3(param3), m_param4(param4) {
   }

   void print_to_log() const;

   std::string format() const;
   /////////////////////////////////////////////
   /// operators
   /////////////////////////////////////////////
	private:
	   /*!
		* \brief operator<<
		*/
		friend std::ostream& operator<<(std::ostream &os, const SynapseInfo& synapse);

		/*!
		 * \brief operator>>
		 */
		friend std::istream& operator>>(std::istream& in, SynapseInfo& synapse);
};

// ///////////////////////////////
// synapse traits				//
// ///////////////////////////////

/**
*  @brief traits for mapping generic synapse info to actual synapse properties
*
*  Always use them for access to generic synapse info if you want to be sure
*  to map values correctly.
*/
template <typename TSynapse = void>
struct synapse_traits
{
	// write access
	static number& loc_coord(SynapseInfo& syn) {return syn.m_locCoords;}
	static unsigned char& type(SynapseInfo& syn) {return syn.m_type;}
	static number& onset(SynapseInfo& syn) {return syn.m_onset;}
	static number& g_max(SynapseInfo& syn) {return syn.m_gMax;}
	static number& v_rev(SynapseInfo& syn) {return syn.m_vRev;}

	// const access
	static number loc_coord(const SynapseInfo& syn) {return syn.m_locCoords;}
	static unsigned char type(const SynapseInfo& syn) {return syn.m_type;}
	static number onset(const SynapseInfo& syn) {return syn.m_onset;}
	static number g_max(const SynapseInfo& syn) {return syn.m_gMax;}
	static number v_rev(const SynapseInfo& syn) {return syn.m_vRev;}
};

template<>
struct synapse_traits<AlphaSynapse> : public synapse_traits<void>
{
	// write access
	static number& tau(SynapseInfo& syn) {return syn.m_tau;}
	static number& freq(SynapseInfo& syn) {return syn.m_param1;}
	static unsigned int& nSpikes(SynapseInfo& syn) {return syn.m_param4;}

	// const access
	static number tau(const SynapseInfo& syn) {return syn.m_tau;}
	static number freq(const SynapseInfo& syn) {return syn.m_param1;}
	static unsigned int nSpikes(const SynapseInfo& syn) {return syn.m_param4;}
	static std::string name() {return std::string("ALPHA_SYNAPSE");}
};

template<>
struct synapse_traits<Exp2Syn> : public synapse_traits<void>
{
	// write access
	static number& tau1(SynapseInfo& syn) {return syn.m_tau;}
	static number& tau2(SynapseInfo& syn) {return syn.m_param3;}
	static number& delay(SynapseInfo& syn) {return syn.m_param1;}
	static number& threshold(SynapseInfo& syn) {return syn.m_param2;}
	static unsigned int& presyn_ind(SynapseInfo& syn) {return syn.m_param4;}

	// const access
	static number tau1(const SynapseInfo& syn) {return syn.m_tau;}
	static number tau2(const SynapseInfo& syn) {return syn.m_param3;}
	static number delay(const SynapseInfo& syn) {return syn.m_param1;}
	static number threshold(const SynapseInfo& syn) {return syn.m_param2;}
	static unsigned int presyn_ind(const SynapseInfo& syn) {return syn.m_param4;}
	static std::string name() {return std::string("BIEXP_SYNAPSE");}

	// functional methods
	static bool activated(const SynapseInfo& syn)
	{
		number on = onset(syn);
		return on == on;
	}
	static void activate(SynapseInfo& syn, number time)
	{
		onset(syn) = time;
	}
	static void deactivate(SynapseInfo& syn)
	{
		onset(syn) = std::numeric_limits<number>::quiet_NaN();
	}
	static number activity_time(const SynapseInfo& syn)
	{
		number t1 = tau1(syn);
		number t2 = tau2(syn);
		number tp = (t1*t2)/(t2-t1) * std::log(t2/t1);
		return tp + 3*t2;
	}
};


} // namespace sh
} // namespace ug
//<! \}

#endif /// __H__UG__SYNAPSE_HANDLER__GRID__SYNAPSE_INFO__
