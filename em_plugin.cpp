/**
 * Plugin of the FE-discretizations of the Maxwell equations.
 * author: D. Logashenko
 */

/* ug headers: */
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

/* general tools for the Nedelec element */
#include "em_material.h"
#include "nedelec_encode.h"
#include "nedelec_gf_user_data.h"
#include "EddyCurrent_E_Nedelec/eddy_current_gf_user_data.h"

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
				.add_method("add", static_cast<void (T::*)(MathVector<dim>&, const char*, const char*)>(&T::add),
							"Sets a constant BC for one component", "Value#Function#Subsets")
				.add_method("add", static_cast<void (T::*)(std::vector<number>&, const char*, const char*)>(&T::add),
							"Sets a constant BC for one component", "Value#Function#Subsets")
				.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >&, const char*, const char*)>(&T::add),
							"Sets a position and time dependent BC for one component", "LuaFunc#Function#Subsets")
#ifdef UG_FOR_LUA
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
			typedef ITransferOperator<TAlgebra> TBase;
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
				.add_method("compute", static_cast<void (T::*)(SmartPtr<GridFunction<TDomain, TAlgebra> >, const char *)>(&T::compute),
							"Evaluates the source field", "GridFunction#Function names")
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
