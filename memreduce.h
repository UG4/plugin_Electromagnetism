/*
 * Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Dmitry Logashenko
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/* ug4 headers */
#include "pcl/pcl.h"
#include "pcl/pcl_communication_structs.h"
#include "common/serialization.h"

/** Implementation of a syncronization of arrays of generic memory structures
 *  distributed among processes on a parallel machine.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__MEM_REDUCE__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__MEM_REDUCE__

namespace ug{
namespace Electromagnetism{

/** Communication policy for the reduction of arrays represented by classes like
 *  std::vector or VariableArray1 using a generic operation. To perform the
 *  the operation, the template parameter class TOp must implement a STATIC
 *  function 'op' with the arguments of type TArray::value_type. During the
 *  reduction, the arguments are initialized with the data on the master and the
 *  current slave.
 *
 * \tparam TArray	type of the array class (e.g. std::vector)
 * \tparam TLayout	type of the parallel layout of the memory
 * \tparam TOp		class of the operation
 */
template <typename TArray, typename TLayout, typename TOp>
class ComPol_MemOp : public pcl::ICommunicationPolicy<TLayout>
{
	typedef ComPol_MemOp<TArray, TLayout, TOp> this_type;
	typedef typename pcl::ICommunicationPolicy<TLayout> base_type;
	
public:
	typedef TArray array_type;
	typedef typename array_type::value_type value_type;
	typedef typename base_type::Layout Layout;
	typedef typename base_type::Interface Interface;

public:
///	Default constructor
	ComPol_MemOp () : m_pMemDst(NULL), m_pMemSrc(NULL) {}

///	Constructor setting the same arrays for the source and the destination
	ComPol_MemOp (TArray* pVec): m_pMemDst(pVec), m_pMemSrc(pVec)	{}

///	Constructor setting the arrays
	ComPol_MemOp (TArray* pMemDst, const TArray* pMemSrc)
		: m_pMemDst(pMemDst), m_pMemSrc(pMemSrc) {}

///	set the source array
	void set_src (const TArray* pMem) {m_pMemSrc = pMem;}

///	set the destination
	void set_dst (TArray* pMemDst) {m_pMemDst = pMemDst;}

/// returns the buffer size
/**
 * This function returns the size of the buffer needed for the communication
 * of passed interface. If the vector has fixed size entries this is just
 * the number of interface entries times the size of the entry. In case
 * of a variable size entry type a negative value is returned to indicate
 * that no buffer size can be determined in advanced.
 *
 * \param[in]	interface	Interface that will communicate
 */
	virtual int get_required_buffer_size (const Interface& interface)
	{
		return interface.size () * sizeof (value_type);
	}

///	writes the interface values into a buffer that will be sent
/**
 * This function collects all entries of the vector into a buffer that
 * are part of the interface.
 *
 * \param[out]		buff		Buffer
 * \param[in]		interface	Interface that will communicate
 */
	virtual bool collect (ug::BinaryBuffer& buff, const Interface& interface)
	{
	//	check that vector has been set
		if (m_pMemSrc == NULL) return false;

	//	loop interface
		for (typename Interface::const_iterator iter = interface.begin ();
				iter != interface.end (); ++iter)
		{
		//	get index
			const size_t index = interface.get_element (iter);

		//	write entry into buffer
			Serialize (buff, (* m_pMemSrc) [index]);
		}
		return true;
	}

///	writes values from a buffer into the interface
/**
 * This function writes the buffer values into the vector.
 *
 * \param[out]		buff		Buffer
 * \param[in]		interface	Interface that communicates
 */
	virtual bool extract (ug::BinaryBuffer& buff, const Interface& interface)
	{
	//	check that vector has been set
		if(m_pMemDst == NULL) return false;

	//	loop interface
		for (typename Interface::const_iterator iter = interface.begin ();
				iter != interface.end (); ++iter)
		{
		//	get index
			const size_t index = interface.get_element (iter);

		//	write entry into the array
			typename TArray::value_type entry;
			Deserialize (buff, entry);
			(* m_pMemDst) [index] = TOp::op ((* m_pMemDst) [index], entry);
		}
		return true;
	}

private:
//	pointer to current vector
	TArray* m_pMemDst;
	const TArray* m_pMemSrc;
};

/** Communication policy for copying arrays represented by classes like
 *  std::vector or VariableArray1.
 *
 * \tparam TArray	type of the array class (e.g. std::vector)
 * \tparam TLayout	type of the parallel layout of the memory
 */
template <typename TArray, typename TLayout>
class ComPol_MemCopy : public pcl::ICommunicationPolicy<TLayout>
{
	typedef ComPol_MemCopy<TArray, TLayout> this_type;
	typedef typename pcl::ICommunicationPolicy<TLayout> base_type;
	
public:
	typedef TArray array_type;
	typedef typename array_type::value_type value_type;
	typedef typename base_type::Layout Layout;
	typedef typename base_type::Interface Interface;

public:
///	Default constructor
	ComPol_MemCopy () : m_pMemDst(NULL), m_pMemSrc(NULL) {}

///	Constructor setting the same arrays for the source and the destination
	ComPol_MemCopy (TArray* pVec): m_pMemDst(pVec), m_pMemSrc(pVec)	{}

///	Constructor setting the arrays
	ComPol_MemCopy (TArray* pMemDst, const TArray* pMemSrc)
		: m_pMemDst(pMemDst), m_pMemSrc(pMemSrc) {}

///	set the source array
	void set_src (const TArray* pMem) {m_pMemSrc = pMem;}

///	set the destination
	void set_dst (TArray* pMemDst) {m_pMemDst = pMemDst;}

/// returns the buffer size
/**
 * This function returns the size of the buffer needed for the communication
 * of passed interface. If the vector has fixed size entries this is just
 * the number of interface entries times the size of the entry. In case
 * of a variable size entry type a negative value is returned to indicate
 * that no buffer size can be determined in advanced.
 *
 * \param[in]	interface	Interface that will communicate
 */
	virtual int get_required_buffer_size (const Interface& interface)
	{
		return interface.size () * sizeof (typename TArray::value_type);
	}

///	writes the interface values into a buffer that will be sent
/**
 * This function collects all entries of the vector into a buffer that
 * are part of the interface.
 *
 * \param[out]		buff		Buffer
 * \param[in]		interface	Interface that will communicate
 */
	virtual bool collect (ug::BinaryBuffer& buff, const Interface& interface)
	{
	//	check that vector has been set
		if (m_pMemSrc == NULL) return false;

	//	loop interface
		for (typename Interface::const_iterator iter = interface.begin ();
				iter != interface.end (); ++iter)
		{
		//	get index
			const size_t index = interface.get_element (iter);

		//	write entry into buffer
			Serialize (buff, (* m_pMemSrc) [index]);
		}
		return true;
	}

///	writes values from a buffer into the interface
/**
 * This function writes the buffer values into the vector.
 *
 * \param[out]		buff		Buffer
 * \param[in]		interface	Interface that communicates
 */
	virtual bool extract (ug::BinaryBuffer& buff, const Interface& interface)
	{
	//	check that vector has been set
		if(m_pMemDst == NULL) return false;

	//	loop interface
		for (typename Interface::const_iterator iter = interface.begin ();
				iter != interface.end (); ++iter)
		{
		//	get index
			const size_t index = interface.get_element (iter);

		//	write entry into array
			Deserialize (buff, (* m_pMemDst) [index]);
		}
		return true;
	}

private:
//	pointer to current vector
	TArray* m_pMemDst;
	const TArray* m_pMemSrc;
};

/** Syncronizes a memory array using a given reduction operation.
 *
 * \tparam TArray	type of the array
 * \tparam TLayout	type of the parallel layout
 * \tparam TOp		class of the reduction operation
 */
template <typename TArray, typename TLayout, typename TOp>
void MemAllReduce
(
	TArray* pMem, ///< the array to syncronize
	const TLayout& masterLayout, ///< the master layout
	const TLayout& slaveLayout, ///< the slave layout
	pcl::InterfaceCommunicator<TLayout>* pCom = NULL ///< (optionally) the communicator
)
{
//	create a new communicator if required
	pcl::InterfaceCommunicator<TLayout> tCom;
	pcl::InterfaceCommunicator<TLayout>& com = (!pCom)? tCom : *pCom;

//	STEP 1: Reduce slave values on the master

	//	create the required communication policies
	ComPol_MemOp<TArray, TLayout, TOp> cpMemOp (pMem);

	//	perform communication
	com.send_data (slaveLayout, cpMemOp);
	com.receive_data (masterLayout, cpMemOp);
	com.communicate ();
	
//	STEP 2: copy master values to slaves

	//	create the required communication policies
	ComPol_MemCopy<TArray, TLayout> cpMemCopy (pMem);

	//	perform communication
	com.send_data (masterLayout, cpMemCopy);
	com.receive_data (slaveLayout, cpMemCopy);
	com.communicate ();
};

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__MEM_REDUCE__

/* End of File */
