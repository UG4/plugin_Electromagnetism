/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/**
 * Plugin of the FE-discretizations of the Maxwell equations.
 */

/* ug headers: */
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

/* general tools for the Nedelec element */
#include "em_material.h"
#include "nedelec_encode.h"
#include "nedelec_gf_user_data.h"
#include "EddyCurrent_E_Nedelec/eddy_current_gf_user_data.h"
#include "EddyCurrent_E_Nedelec/eddy_current_cmd.h"

/* discretizations' headers: */
#include "EddyCurrent_E_Nedelec/eddy_current_e_nedelec.h"

/* hybrid smoother's header: */
#include "EddyCurrent_E_Nedelec/hiptmair_hybrid_smoother.h"

/* transfer operators header: */
#include "nedelec_transfer.h"

/* transfer operators: */
#include "nedelec_dirichlet.h"

/* projection: */
#include "nedelec_project.h"

/* divergence-free sources */
#include "nedelec_source.h"

/* further tools */
#include "nedelec_aux_cmd.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace Electromagnetism{

struct Functionality
{
	/**
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg	registry
	 * @param grp	group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain(Registry& reg, string grp)
	{
		string suffix = GetDomainSuffix<TDomain>();
		string tag = GetDomainTag<TDomain>();
		
	// Parameter specification tools
		{
			typedef EMaterial<TDomain> T;
			string name = string("EMaterial").append(suffix);
			reg.add_class_<T> (name, grp)
				.template add_constructor<void (*) (ConstSmartPtr<TDomain>)>("Domain")
				
				.add_method("add", static_cast<void (T::*) (const char*, number, number)>(&T::add),
							"Adds parameters to a subset", "subset(s)#magn. permeability#electr. conductivity")
				.add_method("add", static_cast<void (T::*) (const char*, number)>(&T::add),
							"Adds parameters of an insulator to a subset", "subset(s)#magn. permeability")
				
				.add_method("close", static_cast<void (T::*) ()>(&T::close),
							"Finalizes the domain description", "")
							
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EMaterial", tag);
		}
	// Further tools
		{
			static const int dim = TDomain::dim;
			string name = string("SubsetIndicatorUserData").append(suffix);
			typedef SubsetIndicatorUserData<TDomain> T;
			typedef UserData<number, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(ConstSmartPtr<TDomain>, const char*)>("Domain#Subsets")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "SubsetIndicatorUserData", tag);
		}
		{
			typedef ug::CPUAlgebra TPotAlgebra;
			typedef ug::GridFunction<TDomain, TPotAlgebra> TPotFct;
			
			reg.add_function("SetSubsetVertVal", static_cast<void (*)(SmartPtr<TPotFct>, const char*, number)>(&SetSubsetVertVal<TPotFct>), grp, "Set a field to a constant on subsets", "GradientGridFunction#Subsets#Value");
		}
	}
	
	/**
	 * Function called for the registration of Domain and Algebra dependent parts.
	 * All Functions and Classes depending on both Domain and Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Domain and Algebra types, based on the current build options.
	 *
	 * @param reg	registry
	 * @param grp	group for sorting of functionality
	 */
	template <typename TDomain, typename TAlgebra>
	static void DomainAlgebra(Registry& reg, string grp)
	{
		static const int dim = TDomain::dim;
		string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
		string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();
		
	// Time-harmonic E-based formulation of the eddy current model:
		{
			typedef EddyCurrent_E_Nedelec<TDomain, TAlgebra> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("EddyCurrent_E_Nedelec").append(suffix);
			reg.add_class_<T, TBase >(name, grp)
				.template add_constructor
					<
						void (*)
						(
							const char*,
							ConstSmartPtr<EMaterial<TDomain> >,
							number
						)
					>("Function(s)#Material data#Frequency")
				.add_method("set_generator_current",
					static_cast<void (T::*)(SmartPtr<GridFunction<TDomain,TAlgebra> >, const char*)>(&T::set_generator_current),
						"Sets the generator current source", "GridFunc#Cmps")
				.add_method("set_generator_current",
					static_cast<void (T::*)(SmartPtr<GridFunction<TDomain,TAlgebra> >, const char*, const char*)>(&T::set_generator_current),
						"Sets the generator current source in subsets", "GridFunc#Cmps#Subssets")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EddyCurrent_E_Nedelec", tag);
		}
	
	// Base class for Dirichlet BC for all type of the discretizations
		{
			typedef EMDirichlet<TDomain, TAlgebra> T;
			typedef IDomainConstraint<TDomain, TAlgebra> TBase;
			string name = string("EMDirichlet").append(suffix);
			reg.add_class_<T, TBase>(name, grp);
			reg.add_class_to_group(name, "EMDirichlet", tag);
		}
	
	// Dirichlet BC for Nedelec-based discretizations
		{
			typedef NedelecDirichletBC<TDomain, TAlgebra> T;
			typedef EMDirichlet<TDomain, TAlgebra> TBase;
			string name = string("NedelecDirichletBC").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*) (const char*)>("Functions")
				.add_method("add_0", static_cast<void (T::*)(const char*)>(&T::add_0),
							"Sets a zero BC on subsets", "Value#Function#Subsets")
#ifdef UG_FOR_LUA
				.add_method("add", static_cast<void (T::*)(MathVector<dim>&, const char*, const char*)>(&T::add),
							"Sets a constant BC for one component", "Value#Function#Subsets")
				.add_method("add", static_cast<void (T::*)(std::vector<number>&, const char*, const char*)>(&T::add),
							"Sets a constant BC for one component", "Value#Function#Subsets")
				.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >&, const char*, const char*)>(&T::add),
							"Sets a position and time dependent BC for one component", "UserData#Function#Subsets")
				.add_method("add", static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add),
							"Sets a position and time dependent BC for one component", "LuaFunc#Function#Subsets")
#endif
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecDirichletBC", tag);
		}
	
	// Hybrid smoother by Hiptmair
		{
			typedef TimeHarmonicNedelecHybridSmoother<TDomain, TAlgebra> T;
			typedef ILinearIterator<typename T::vector_type> TBase;
			string name = string("HiptmairHybridSmoother").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor
					<
						void (*)
						(
							SmartPtr<ApproximationSpace<TDomain> >,
							SmartPtr<ILinearIterator<typename T::vector_type> >,
							SmartPtr<ILinearIterator<typename T::pot_vector_type> >
						)
					>("VertexApproxSpace#EdgeSmoother#VertexSmoother")
				.add_method("set_Dirichlet", static_cast<void (T::*)(SmartPtr<EMDirichlet<TDomain, TAlgebra> >)>(&T::set_Dirichlet),
							"Sets the object of the Dirichlet BC", "Dirichlet BC")
				.add_method("set_skip_edge_smoother", static_cast<void (T::*)(bool)>(&T::set_skip_edge_smoother),
							"Whether to skip the edge smoother (for debugging)", "Flag")
				.add_method("set_skip_vertex_smoother", static_cast<void (T::*)(bool)>(&T::set_skip_vertex_smoother),
							"Whether to skip the vertex smoother (for debugging)", "Flag")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "HiptmairHybridSmoother", tag);
		}
	
	// Transfer operators for the Whitney-1 elements
		{
			typedef NedelecTransfer<TDomain, TAlgebra> T;
			typedef ITransferOperator<TDomain, TAlgebra> TBase;
			string name = string("NedelecTransfer").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor
					<
						void (*)
						(
							SmartPtr<ApproximationSpace<TDomain> >
						)
					>("ApproxSpace")
				.add_method("add_constraint", &T::add_constraint)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecTransfer", tag);
		}
		
	// Projection solution to the divergence-free space
		{
			typedef NedelecProject<TDomain, TAlgebra> T;
			string name = string("NedelecProject").append(suffix);
			reg.add_class_<T>(name, grp)
				.template add_constructor
					<
						void (*)
						(
							SmartPtr<EMaterial<TDomain> >,
							SmartPtr<ApproximationSpace<TDomain> >,
							SmartPtr<ILinearOperatorInverse<typename NedelecProject<TDomain, TAlgebra>::pot_vector_type> >
						)
					>("Matherial data#Vert.-based ApproxSpace#Vert.-based LinSolver")
				.add_method("set_Dirichlet", static_cast<void (T::*)(SmartPtr<EMDirichlet<TDomain, TAlgebra> >)>(&T::set_Dirichlet),
							"Sets the object of the Dirichlet BC", "Dirichlet BC")
				.add_method("apply", static_cast<void (T::*)(SmartPtr<GridFunction<TDomain, TAlgebra> >, const char *)>(&T::apply),
							"Projects given functions", "GridFunction#Function names")
				.add_method("compute_div", static_cast<void (T::*)(SmartPtr<GridFunction<TDomain, TAlgebra> >, const char *,SmartPtr<GridFunction<TDomain, typename T::TPotAlgebra> >)>(&T::compute_div),
							"Compute weak div in insulators", "GridFunction for u#Function names#GridFunction for div")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecProject", tag);
		}
	
	// Computation of the Whitney-1 DoFs for a given function
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			typedef ug::CPUAlgebra TPotAlgebra;
			typedef ug::GridFunction<TDomain, TPotAlgebra> TPotFct;
			
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*, const char*, number)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "Data#GridFunction#Component#Subsets#Time");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*, number)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "Data#GridFunction#Component#Time");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*, const char*)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "Data#GridFunction#Component#Subsets");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "Data#GridFunction#Component");
			
			#ifdef UG_FOR_LUA
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*, number)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "LuaFunction#GridFunction#Component#Subsets#Time");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, number)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "LuaFunction#GridFunction#Component#Time");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "LuaFunction#GridFunction#Component#Subsets");
			reg.add_function("ComputeNedelecDoFs", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*)>(&ComputeNedelecDoFs<TFct>), grp, "Nedelec DoFs for given vector field", "LuaFunction#GridFunction#Component");
			#endif
			
			reg.add_function("NedelecGradPotential", static_cast<void (*)(SmartPtr<TPotFct>, SmartPtr<TFct>, const char*)>(&NedelecGradPotential<TPotFct, TFct>), grp, "Nedelec DoFs for the gradient of a given field", "GradientGridFunction#ResultGridFunction#Component");
		}
	
	//	Computation of divergence-free sources
		{
			typedef NedelecLoopCurrent<TDomain, TAlgebra> T;
			string name = string("NedelecLoopCurrent").append(suffix);
			reg.add_class_<T>(name, grp)
				.template add_constructor
					<
						void (*)
						(
							const char *, const char *, const char *,
							SmartPtr<ApproximationSpace<TDomain> >,
							SmartPtr<ILinearOperatorInverse<typename NedelecProject<TDomain, TAlgebra>::pot_vector_type> >
						)
					>("Source subsets#Pos. dir. subsets#Cut subsets#Vert. approx. space#Lin. solver for potential")
				.add_method("set", static_cast<void (T::*)(const char *, number)>(&T::set),
							"Sets the electric current", "Component#Value")
				.add_method("compute", static_cast<void (T::*)(SmartPtr<GridFunction<TDomain, TAlgebra> >)>(&T::compute),
							"Evaluates the source field", "GridFunction")
				.add_method("subsets", static_cast<std::string (T::*)()>(&T::subsets),
							"Returns the source's subsets", "")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecLoopCurrent", tag);
		}
	
	//	Computation of the vector and curl fields for a given Nedelec-element based grid function, etc
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("NedelecGridFunctionData").append(suffix);
			typedef NedelecGridFunctionData<TFct> T;
			typedef UserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Components")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecGridFunctionData", tag);
		}
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("NedelecCurlData").append(suffix);
			typedef NedelecCurlData<TFct> T;
			typedef UserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Components")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecCurlData", tag);
		}
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("NedelecSigmaEData").append(suffix);
			typedef NedelecSigmaEData<TFct> T;
			typedef UserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*, SmartPtr<EMaterial<TDomain> >)>("GridFunction#Components#Materials")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NedelecSigmaEData", tag);
		}
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("EddyCurrentHeat").append(suffix);
			typedef EddyCurrentHeat<TFct> T;
			typedef UserData<number, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*, SmartPtr<EMaterial<TDomain> >)>("GridFunction#Components#Materials")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EddyCurrentHeat", tag);
		}
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("EddyCurrentReBofEUserData").append(suffix);
			typedef EddyCurrentReBofEUserData<TFct> T;
			typedef UserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*, number)>("GridFunction#Components#Frequency")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EddyCurrentReBofEUserData", tag);
		}
		{
			static const int dim = TDomain::dim;
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			string name = string("EddyCurrentImBofEUserData").append(suffix);
			typedef EddyCurrentImBofEUserData<TFct> T;
			typedef UserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<TFct>, const char*, number)>("GridFunction#Components#Frequency")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EddyCurrentImBofEUserData", tag);
		}
		
	// Computation of various values
		{
			typedef ug::GridFunction<TDomain, TAlgebra> TFct;
			
			reg.add_function
			(
				"CalcPower",
				static_cast
				<
					void (*)
					(
						SmartPtr<TFct> spJGGF,
						const char* JG_cmps,
						const char* JG_ss,
						SmartPtr<TFct> spEGF,
						const char* E_cmps
					)
				>
				(&CalcPower<TFct>),
				grp,
				"Power of the electromagnetic field (up to the contribution of the boundary)",
				"GeneratorCurrent#cmps#SubSets#ElectricField#cmps"
			);
			
			reg.add_function
			(
				"CalcEMF",
				static_cast
				<
					void (*)
					(
						SmartPtr<TFct>,
						const char*,
						const char*,
						const std::vector<number>&,
						const std::vector<number>&,
						const size_t,
						const std::vector<number>&
					)
				>
				(&CalcEMF<TFct>),
				grp,
				"Magnetic flux through a cylindric coil",
				"ElectricField#cmps#subsets#normal#basePnt#numWindings#windingSize"
			);
			
			reg.add_function
			(
				"ComputeFlux",
				static_cast
				<
					void (*)
					(
						SmartPtr<TFct>,
						const char*,
						const char*,
						const char*
					)
				>
				(&ComputeFlux<TFct>),
				grp,
				"Flux of a Nedelec vector field through a surface",
				"Field#cmp#vol.subsets#surface subsets"
			);
		}
	};
	
};

} // end namespace Electromagnetism

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_Electromagnetism(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/Electromagnetism");
	typedef Electromagnetism::Functionality Functionality;

	try{
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace ug

/* End of File */
